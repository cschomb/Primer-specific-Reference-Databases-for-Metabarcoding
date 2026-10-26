import argparse
import gzip
import re
import logging
import pandas as pd
from tqdm import tqdm

RANK_MAP = {
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species",
}
TARGET_RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return open(path, mode="rt", encoding="utf-8", errors="replace")

def parse_taxonomy_from_header(header: str):
    m = re.search(r"tax=([^\s]+)", header)
    if not m:
        return None
    tax_str = m.group(1)
    parts = [p for p in tax_str.split(",") if p]
    out = {}
    for part in parts:
        if ":" not in part:
            continue
        code, value = part.split(":", 1)
        code = code.strip().lower()
        value = value.strip()
        rank_full = RANK_MAP.get(code)
        if rank_full and value:
            out[rank_full] = value
    if not any(r in out for r in TARGET_RANKS):
        return None
    return out

def iter_headers(path: str):
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                yield line[1:]
            elif "tax=" in line:
                # fallback: non-standard header with tax=
                yield line

def records_from_fasta_with_progress(path: str):
    logging.info(f"Loading records from file: {path}.\n ")
    headers = list(iter_headers(path))  # preload all headers
    records = []
    for header in tqdm(headers, desc="Processing records", unit="records"):
        tax = parse_taxonomy_from_header(header)
        if tax is None:
            continue
        rec = {rank: tax.get(rank) for rank in TARGET_RANKS}
        records.append(rec)
    print()
    return records

def autosize_columns(worksheet, df, max_width: int = 50):
    for i, col in enumerate(df.columns):
        series = df[col].astype(str)
        max_len = max([len(str(col))] + [len(str(x)) for x in series])
        worksheet.set_column(i, i, min(max_len + 2, max_width))

def write_excel(df: pd.DataFrame, out_path: str):
    if df.empty:
        raise ValueError("No taxonomy records found in the input. Check header format (tax=...).")
    logging.info(f"Writing processed records to file.\n ")
    df = df.copy()
    df[TARGET_RANKS] = df[TARGET_RANKS].fillna("Unassigned")
    with pd.ExcelWriter(out_path) as writer:
        for i, rank in enumerate(TARGET_RANKS):
            cols = TARGET_RANKS[: i + 1]
            grouped = (
                df.groupby(cols, dropna=False)
                  .size()
                  .reset_index(name="Count")
                  .sort_values(["Count"] + cols, ascending=[False] + [True] * len(cols), kind="mergesort")
            )
            grouped.to_excel(writer, sheet_name=rank, index=False)

            ws = writer.sheets[rank]
            try:
                ws.freeze_panes(1, 0)
                autosize_columns(ws, grouped)
            except Exception:
                pass
    logging.info(f"Wrote Excel workbook to: {out_path}")

def main():
    ap = argparse.ArgumentParser(description="Summarize taxonomy from FASTA headers into an Excel workbook. Format of the FASTA headers needs to be in the following format: >ACCESSION;tax=p:Phylum,c:CLASS,o:ORDER,f:FAMILY,g:GENUS,s:SPECIES")
    ap.add_argument("fasta", help="Input FASTA file (can be .gz)")
    ap.add_argument("outxlsx", help="Output Excel file path (.xlsx)")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    records = records_from_fasta_with_progress(args.fasta)
    if not records:
        raise SystemExit("No headers with 'tax=' found or unable to parse taxonomy. Please check input format.")
    df = pd.DataFrame.from_records(records, columns=TARGET_RANKS)
    write_excel(df, args.outxlsx)

if __name__ == "__main__":
    main()

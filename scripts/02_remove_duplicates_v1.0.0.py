import argparse
import re
from collections import defaultdict
from tqdm import tqdm

# Match entries like "ACCESSION;tax=...<taxonomy...>" up to the next ", ACCESSION;tax=" or end of line
ENTRY_RE = re.compile(r'([A-Za-z0-9_.]+);tax=(.*?)(?=(?:,\s*[A-Za-z0-9_.]+;tax=)|$)')

def normalize_taxonomy(taxstr: str) -> str:
    """
    Normalize the taxonomy substring (everything after ';tax=') so that
    trivial formatting differences don't prevent matches.
    Strategy:
      - strip leading/trailing whitespace
      - drop any trailing commas
      - remove all internal whitespace
    """
    s = taxstr.strip()
    # Remove trailing commas and spaces
    while s.endswith(",") or (s and s[-1].isspace()):
        s = s.rstrip().rstrip(",")
    # Remove all whitespace to avoid comma/space differences affecting comparison
    s = re.sub(r"\s+", "", s)
    return s

def iter_entries(rest: str):
    """
    Yield (accession, normalized_taxonomy, raw_taxonomy) for the comma-separated portion after the leading count.
    """
    for m in ENTRY_RE.finditer(rest):
        acc = m.group(1).strip()
        tax_raw = m.group(2)
        tax_key = normalize_taxonomy(tax_raw)
        yield acc, tax_key, tax_raw.strip()

def process_file(infile: str, prefer_first: bool = True):
    accepted_set = set()
    accepted_order = []  # preserve acceptance order

    remove_set = set()
    remove_order = []    # preserve first-removal order

    actions = []         # per-line report

    with open(infile, "r", encoding="utf-8") as fh:
        total_lines = sum(1 for _ in fh)

    with open(infile, "r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(tqdm(fh, total=total_lines, desc="Processing file lines"), start=1):
            line = raw.strip()
            if not line:
                continue

            # Split off leading count (tab or whitespace)
            try:
                n_str, rest = line.split("\t", 1)
            except ValueError:
                parts = line.split(None, 1)
                if len(parts) != 2:
                    continue
                n_str, rest = parts

            # Parse entries
            entries = list(iter_entries(rest))
            if not entries:
                continue

            # Group by full taxonomy key (already normalized)
            groups = defaultdict(list)   # tax_key -> list of accessions (input order)
            tax_readable = {}            # tax_key -> original (trimmed) taxonomy for reporting
            for acc, tax_key, tax_raw in entries:
                groups[tax_key].append(acc)
                if tax_key not in tax_readable:
                    tax_readable[tax_key] = tax_raw

            # Process each taxonomy group in this line
            for tax_key, accs in groups.items():
                notes = []
                preaccepted = [a for a in accs if a in accepted_set]

                if preaccepted:
                    kept = list(preaccepted)
                    removed_here = []
                    for a in accs:
                        if a in accepted_set:
                            continue
                        if a not in remove_set:
                            remove_set.add(a)
                            remove_order.append(a)
                            removed_here.append(a)
                    if len(preaccepted) > 1:
                        notes.append("multiple_preaccepted_in_group")

                    actions.append({
                        "line_no": line_no,
                        "taxonomy": tax_readable[tax_key],
                        "kept": ",".join(kept) if kept else "",
                        "removed": ",".join(removed_here) if removed_here else "",
                        "note": "|".join(notes) if notes else ""
                    })
                else:
                    keeper = accs[0] if prefer_first else sorted(accs)[0]
                    kept = [keeper]
                    removed_here = []

                    if keeper in remove_set:
                        remove_set.remove(keeper)
                        try:
                            remove_order.remove(keeper)
                        except ValueError:
                            pass
                        notes.append("rescued_previously_marked")

                    if keeper not in accepted_set:
                        accepted_set.add(keeper)
                        accepted_order.append(keeper)

                    for a in accs[1:]:
                        if a in accepted_set:
                            notes.append("unexpected_preaccepted_duplicate")
                            continue
                        if a not in remove_set:
                            remove_set.add(a)
                            remove_order.append(a)
                            removed_here.append(a)

                    actions.append({
                        "line_no": line_no,
                        "taxonomy": tax_readable[tax_key],
                        "kept": ",".join(kept),
                        "removed": ",".join(removed_here) if removed_here else "",
                        "note": "|".join(sorted(set(notes))) if notes else ""
                    })

    return accepted_order, remove_order, actions

def main():
    ap = argparse.ArgumentParser(description="Keep one accession per unique taxonomy per line; never remove previously accepted accessions.")
    ap.add_argument("input", help="Input file")
    ap.add_argument("--out-prefix", default="out", help="Output prefix (default: out)")
    ap.add_argument("--prefer", choices=["first","lex"], default="first",
                    help="Keeper selection within a line's taxonomy group: first (input order) or lex (lexicographic). Default: first")
    args = ap.parse_args()

    prefer_first = (args.prefer == "first")
    accepted, removed, actions = process_file(args.input, prefer_first=prefer_first)

    with open(f"{args.out_prefix}.accepted.txt", "w", encoding="utf-8") as fa:
        for acc in tqdm(accepted, desc="Writing accepted"):
            fa.write(acc + "\n")

    with open(f"{args.out_prefix}.remove.txt", "w", encoding="utf-8") as fr:
        for acc in tqdm(removed, desc="Writing removed"):
            fr.write(acc + "\n")

    with open(f"{args.out_prefix}.per_line_actions.tsv", "w", encoding="utf-8") as fp:
        fp.write("line_no\ttaxonomy\tkept\tremoved\tnote\n")
        for row in tqdm(actions, desc="Writing per-line actions"):
            fp.write(f"{row['line_no']}\t{row['taxonomy']}\t{row['kept']}\t{row['removed']}\t{row['note']}\n")

    print(f"Wrote {len(accepted)} accepted to {args.out_prefix}.accepted.txt")
    print(f"Wrote {len(removed)} to-remove to {args.out_prefix}.remove.txt")
    print(f"Wrote per-line details to {args.out_prefix}.per_line_actions.tsv")

if __name__ == "__main__":
    main()

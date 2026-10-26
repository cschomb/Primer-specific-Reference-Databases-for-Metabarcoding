import os
import argparse
import logging
import pandas as pd
import json
from tqdm import tqdm
tqdm.pandas()

#-- ASSIGN CATEGORY, TAKING INTO ACCOUNT THE TOTAL NUMBER OF SPECIES SEQUENCES IN THE BLAST DURING THE categorize_and_write STEP --#
def assign_category(row):
    count = row['species_count']
    levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    labels = ['P', 'C', 'O', 'F', 'G', 'S']
    label_to_level = dict(zip(labels, levels))

    highest_diff_level = row['highest_diff_level']
    diff_json = row.get('difference_values')
    q_species = str(row['q_species'])
    q_level_value = None

    if not diff_json or diff_json == '{}' or pd.isna(diff_json):
        n_main = 0
        total = 0
    else:
        try:
            diff_dict = json.loads(diff_json)
        except Exception:
            return 'unknown'
        total = sum(diff_dict.values())

        if not highest_diff_level or pd.isna(highest_diff_level):
            n_main = 0
            q_level_value = None
        else:
            col_name = f"q_{label_to_level.get(highest_diff_level, 'species')}"
            q_level_value = str(row.get(col_name, ''))
            n_main = diff_dict.get(q_level_value, 0)

    if count < 3:
        if pd.isna(highest_diff_level) or highest_diff_level == '':
            return 'singletons'
        else:
            return 'low_coverage'
    else:
        if not highest_diff_level or pd.isna(highest_diff_level):
            return 'accepted'
        if total == 0:
            return 'check'
        frac_main = n_main / total if total > 0 else 0

        if frac_main >= 0.75:
            return 'accepted'
        elif frac_main <= 0.25:
            if highest_diff_level == 'S':
                return 'check'
            else:
                return 'refused'
        else:
            return 'check'

#-- WRITE CATEGORISED SEQUENCES TO FILES --#
def categorize_and_write(summary_csv, species_count_file, output_path, out_prefix='category_'):
    summary_df = pd.read_csv(summary_csv)
    scount = pd.read_csv(species_count_file, sep='\t', names=['q_species', 'species_count'])

    merged = pd.merge(summary_df, scount, on='q_species', how='left')

    categories = []
    for _, row in tqdm(merged.iterrows(), total=len(merged), desc="Assigning categories"):
        categories.append(assign_category(row))
    merged['category'] = categories

    for category, group in merged.groupby('category'):
        file_name = f"{out_prefix}{category}.csv"
        out_file = os.path.join(output_path, file_name)
        group.to_csv(out_file, index=False)
        logging.info(f"Written {len(group)} entries to {out_file}")

#-- COMPARES LEVELS OF TAXONOMIC DIFFERENCE FOR THE BLAST RESULTS DURING THE analyze_per_qacc STEP --#
def compare_taxonomy(row):
    levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    labels = ['P', 'C', 'O', 'F', 'G', 'S']

    diffs = []
    for lvl, label in zip(levels, labels):
        q_val = row[f'q_{lvl}']
        s_val = row[f's_{lvl}']
        diff = 0 if q_val == s_val else 1
        diffs.append(f"{label}:{diff}")
    return ", ".join(diffs)

#-- ANALYSES THE BLAST RESULTS PER QUERY AND WRITES SUMMARY TABLE --#
def analyze_per_qacc(csv_file, output_path):
    filtered_df = pd.read_csv(csv_file)
    summary_rows = []

    levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    labels = ['P', 'C', 'O', 'F', 'G', 'S']
    label_to_rank = dict(zip(labels, levels))

    for q_acc, group in filtered_df.groupby('q_acc'):
        q_tax = group.iloc[0][[f'q_{lvl}' for lvl in levels]].to_dict()
        group_size = len(group)

        diff_strings = group.apply(compare_taxonomy, axis=1)
        score_interval = group["score_interval"].iloc[0]

        diff_counts = dict.fromkeys(labels, 0)
        for diff_str in diff_strings:
            parts = diff_str.split(", ")
            for part in parts:
                label, val = part.split(":")
                diff_counts[label] += int(val)

        diff_counts_str = ", ".join([f"{label}:{diff_counts[label]}" for label in labels])

        highest_diff_level = None
        for label in labels:
            if diff_counts[label] > 0:
                highest_diff_level = label
                break

        if highest_diff_level is not None:
            diff_rank = label_to_rank[highest_diff_level]
            subject_values = group[f's_{diff_rank}'].dropna().astype(str).tolist()
            all_values = subject_values + [q_tax[f'q_{diff_rank}']]

            from collections import Counter
            counts = dict(Counter(all_values))
            difference_values = counts
        else:
            highest_diff_level = ''
            difference_values = {}

        difference_values_str = json.dumps(difference_values)

        summary_rows.append({
            'q_acc': q_acc,
            **{f'q_{lvl}': q_tax[f'q_{lvl}'] for lvl in levels},
            'diff_counts': diff_counts_str,
            'highest_diff_level': highest_diff_level,
            'group_size': group_size,
            'score_interval': score_interval,
            'difference_values': difference_values_str
        })

    summary_df = pd.DataFrame(summary_rows)
    output_file = os.path.join(output_path, 'summarised_blast_results.csv')
    summary_df.to_csv(output_file, index=False)
    logging.info(f"Summary data per q_acc saved to {output_file}")

#-- EXTRACTS THE TAXONOMY OF THE SEQUENCES DURING THE process_and_log_group STEP --#
def extract_taxonomy(tax_str):
    if pd.isnull(tax_str):
        return [None] * 6
    tax_parts = tax_str.replace("tax=", "").split(',')
    return [
        tax_parts[i].split(':')[1] if len(tax_parts) > i else None
        for i in range(6)
    ]

#-- FILTERS BLAST INPUT FILE IN CHUNKS --#
def filter_blast(input_file, output_path, length_threshold, chunksize=10000):
    logging.info(f"Loading BLAST file: {input_file}.")
    output_csv = os.path.join(output_path, 'filtered_blast_results.csv')
    species_count_file = os.path.join(output_path, 'species_count.txt')

    species_count = {}

    if os.path.exists(output_csv):
        os.remove(output_csv)

    curr_q_acc = None
    buffer = []
    column_order = [
        'q_acc', 'q_phylum', 'q_class', 'q_order', 'q_family',
        'q_genus', 'q_species', 's_acc', 's_phylum', 's_class',
        's_order', 's_family', 's_genus', 's_species', 'pident',
        'length', 'mismatches', 'gapopens', 'qstart', 'qend',
        'sstart', 'send', 'evalue', 'bitscore', 'qcov', 'score', 'score_interval'
    ]
    write_header = True

    import math
    try:
        total_lines = sum(1 for _ in open(input_file))
        n_chunks = math.ceil(total_lines / chunksize)
    except Exception:
        n_chunks = None

    chunk_number = 0
    logging.info("Starting to process input file in chunks.")

    def process_and_log_group(group_df, q_acc_name):
        n_input = len(group_df)
        group_df = group_df[group_df['length'] >= length_threshold].copy()
        if group_df.empty:
            tqdm.write(f"Processing q_acc: {q_acc_name} ({n_input} hits), kept 0 after filtering.")
            return None
        group_df['pident'] = group_df['pident'].round(3)
        group_df['score'] = (group_df['pident'] + (group_df['qcov'] / 5)).round(3)
        group_df = group_df.sort_values('score', ascending=False)

        max_score = group_df['score'].max()
        min_score = max_score - 3
        selected_group = group_df[group_df['score'] >= min_score].copy()

        n_kept = len(selected_group)
        if n_kept > 0:
            interval_str = f"{max_score}-{round(min_score, 3)}"
            selected_group['score_interval'] = interval_str
            tqdm.write(f"Processing q_acc: {q_acc_name} ({n_input} hits), kept {n_kept} after filtering, score interval {interval_str}")
        else:
            tqdm.write(f"Processing q_acc: {q_acc_name} ({n_input} hits), kept 0 after filtering.")
        return selected_group

    for chunk in tqdm(pd.read_csv(input_file, sep='\t', header=None, chunksize=chunksize), total=n_chunks, desc="Filtering BLAST chunks"):
        chunk_number += 1
        tqdm.write(f"Processing chunk #{chunk_number}")

        chunk.columns = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatches', 'gapopens',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcov'
        ]
        chunk[['q_acc', 'q_tax']] = chunk['qseqid'].str.split(';', n=1, expand=True)
        chunk[['s_acc', 's_tax']] = chunk['sseqid'].str.split(';', n=1, expand=True)
        chunk[['q_phylum', 'q_class', 'q_order', 'q_family', 'q_genus', 'q_species']] = chunk['q_tax'].apply(extract_taxonomy).apply(pd.Series)
        chunk[['s_phylum', 's_class', 's_order', 's_family', 's_genus', 's_species']] = chunk['s_tax'].apply(extract_taxonomy).apply(pd.Series)
        chunk.drop(['qseqid', 'sseqid', 'q_tax', 's_tax'], axis=1, inplace=True)
        chunk['length'] = pd.to_numeric(chunk['length'], errors='coerce')

        for _, row in chunk.iterrows():
            this_q_acc = row['q_acc']
            if curr_q_acc is None:
                curr_q_acc = this_q_acc
            if this_q_acc != curr_q_acc:
                filtered = process_and_log_group(pd.DataFrame(buffer), curr_q_acc)
                if filtered is not None and not filtered.empty:
                    filtered = filtered[column_order]
                    counts = filtered.groupby('q_species')['q_acc'].nunique()
                    for sp, count in counts.items():
                        species_count[sp] = species_count.get(sp, 0) + count
                    filtered.to_csv(output_csv, mode='a', header=write_header, index=False)
                    if write_header:
                        write_header = False
                buffer = []
                curr_q_acc = this_q_acc
            buffer.append(row)

    if buffer:
        filtered = process_and_log_group(pd.DataFrame(buffer), curr_q_acc)
        if filtered is not None and not filtered.empty:
            filtered = filtered[column_order]
            counts = filtered.groupby('q_species')['q_acc'].nunique()
            for sp, count in counts.items():
                species_count[sp] = species_count.get(sp, 0) + count
            filtered.to_csv(output_csv, mode='a', header=write_header, index=False)
            write_header = False

    with open(species_count_file, 'w') as f:
        for sp, ct in species_count.items():
            f.write(f"{sp}\t{ct}\n")

def main():
    parser = argparse.ArgumentParser(description='Analyze BLAST output file in outfmt 6 format.')

    # --- Input mode: either raw BLAST (-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs") or pre-filtered files ---
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        'input_file',
        nargs='?',
        help='Path to the raw BLAST output file in outfmt 6 format (-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs").'
    )
    input_group.add_argument(
        '--resume',
        nargs=2,
        metavar=('FILTERED_BLAST_CSV', 'SPECIES_COUNT_TXT'),
        help=(
            'Skip filtering and start from pre-existing files. '
            'Provide paths to filtered_blast_results.csv and species_count.txt.'
        )
    )

    parser.add_argument(
        '-o', '--output_dir',
        default='.',
        help='Output directory for results (default: current directory).'
    )
    parser.add_argument(
        '-t', '--length_threshold',
        type=int,
        default=200,
        help='Minimum alignment length to filter results (default: 200). Only used without --resume.'
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    if args.resume:
        # --- Resume mode: skip filtering ---
        filtered_csv_file, species_count_file = args.resume

        if not os.path.exists(filtered_csv_file):
            logging.error(f"Filtered BLAST CSV not found: {filtered_csv_file}")
            raise SystemExit(1)
        if not os.path.exists(species_count_file):
            logging.error(f"Species count file not found: {species_count_file}")
            raise SystemExit(1)

        logging.info(f"Resume mode: using {filtered_csv_file} and {species_count_file}")

    else:
        # --- Full pipeline: start from raw BLAST file ---
        if not os.path.exists(args.input_file):
            logging.error(f"Input BLAST file not found: {args.input_file}")
            raise SystemExit(1)

        filter_blast(args.input_file, output_dir, args.length_threshold)

        filtered_csv_file = os.path.join(output_dir, 'filtered_blast_results.csv')
        species_count_file = os.path.join(output_dir, 'species_count.txt')

    # --- Steps shared by both modes ---
    analyze_per_qacc(filtered_csv_file, output_dir)

    summarised_csv_file = os.path.join(output_dir, 'summarised_blast_results.csv')
    categorize_and_write(summarised_csv_file, species_count_file, output_dir, out_prefix='blast_')

    logging.info("FINISHED!")


if __name__ == '__main__':
    main()

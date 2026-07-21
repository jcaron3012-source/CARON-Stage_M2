import sys
import csv
from pathlib import Path
from collections import defaultdict

'''For each FASTA file in a repository, calculates the percentage of identity of each species compared with the mouse sequence.
Returns a matrix species x sequences with % id.''''

def parse_fasta(filepath):
    sequences = {}
    order = []
    current_name = None
    current_seq = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
                order.append(current_name)
            elif line:
                current_seq.append(line.strip())
    if current_name is not None:
        sequences[current_name] = "".join(current_seq)
    return sequences, order

def check_alignment(sequences):
    '''Check if every sequence has the same number of bases'''
    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) > 1:
        sys.exit(1)
    return lengths.pop()

def compute_identity(ref_seq, query_seq):
    '''Returns a % of identity for each sequence'''
    compared = 0
    identical = 0
    for r, q in zip(ref_seq, query_seq):
        if r == "-":
            continue
        compared += 1
        if r.upper() == q.upper():
            identical += 1
    return (identical / compared * 100) if compared else 0.0

def process_file(fasta_path, ref_name=None):
    sequences, order = parse_fasta(fasta_path)
    check_alignment(sequences)
    if ref_name is  None:
        sys.exit(1)
    if ref_name not in sequences:
        sys.exit(1)
    ref_seq = sequences[ref_name]
    results = {}
    for name in order:
        if name == ref_name:
            continue
        results[name] = compute_identity(ref_seq, sequences[name])
    mean = sum(results.values()) / len(results) if results else 0.0
    return results, mean

def main():
    if len(sys.argv) < 2:
        sys.exit(1)

    fasta_files = [Path(p) for p in sys.argv[1:]]
    all_species = set()
    per_file_results = {}
    per_file_mean = {}

    for fasta in fasta_files:
        print(f"Processing {fasta}")
        results, mean = process_file(fasta, "Mus_musculus")
        per_file_results[fasta.stem] = results
        per_file_mean[fasta.stem] = mean
        all_species.update(results.keys())

    all_species = sorted(all_species)
    out_csv = Path("identity_matrix.csv")

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["fasta"] + all_species + ["MEAN"]
        writer.writerow(header)

        species_values = defaultdict(list)

        for fasta_name in per_file_results:
            row = [fasta_name]
            results = per_file_results[fasta_name]
            for sp in all_species:
                if sp in results:
                    val = round(results[sp], 4)
                    row.append(val)
                    species_values[sp].append(val)
                else:
                    row.append("NA")
            row.append(round(per_file_mean[fasta_name], 4))
            writer.writerow(row)

        mean_row = ["MEAN"]
        for sp in all_species:
            vals = species_values[sp]
            if vals:
                mean_row.append(round(sum(vals) / len(vals), 4))
            else:
                mean_row.append("NA")
        mean_row.append("")
        writer.writerow(mean_row)

if __name__ == "__main__":
    main()

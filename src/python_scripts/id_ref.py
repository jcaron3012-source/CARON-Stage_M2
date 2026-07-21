import sys
import csv
from pathlib import Path

"""Computes the % identity of each sequence vs a reference."""

def parse_fasta(filepath):
    """Parse a multi-FASTA file. Returns a dict {name: sequence} and an ordered list of names."""
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
                current_name = line[1:].split()[0]  # First word after ">"
                current_seq = []
                order.append(current_name)
            elif line:
                current_seq.append(line)

    if current_name is not None:
        sequences[current_name] = "".join(current_seq)

    return sequences, order


def check_alignment(sequences):
    """Check that all sequences have the same length (i.e. are aligned)."""
    lengths = {name: len(seq) for name, seq in sequences.items()}
    unique_lengths = set(lengths.values())
    if len(unique_lengths) > 1:
        print("WARNING: sequences do not all have the same length!")
        for name, length in lengths.items():
            print(f"  {name}: {length} bp")
        sys.exit(1)
    return list(unique_lengths)[0]


def compute_identity(ref_seq, query_seq):
    """
    Compute the % identity between ref_seq and query_seq.
    Positions where the reference has a gap are ignored.
    Returns (% identity, positions compared, identical positions).
    """
    compared = 0
    identical = 0

    for ref_base, query_base in zip(ref_seq, query_seq):
        if ref_base == "-":
            continue  # Skip gap positions in the reference
        compared += 1
        if ref_base.upper() == query_base.upper():
            identical += 1

    pct = (identical / compared * 100) if compared > 0 else 0.0
    return pct, compared, identical


def main():
    #Arguments
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    fasta_path = Path(sys.argv[1])
    if not fasta_path.exists():
        print(f"Error: file not found: {fasta_path}")
        sys.exit(1)

    #Parsing
    print(f"\nReading: {fasta_path.name}")
    sequences, order = parse_fasta(fasta_path)
    print(f"  {len(sequences)} sequences found.\n")

    #List available sequences
    print("Available sequences:")
    for i, name in enumerate(order, 1):
        print(f"  [{i:2d}] {name}")

    #Choose reference
    if len(sys.argv) >= 3:
        ref_name = sys.argv[2]
        if ref_name not in sequences:
            print(f"\nError: reference '{ref_name}' not found in the file.")
            sys.exit(1)
        print(f"\nReference (from argument): {ref_name}")
    else:
        print("\nEnter the exact name of the reference sequence")
        print("(or its number from the list above): ", end="")
        user_input = input().strip()

        if user_input.isdigit():
            idx = int(user_input) - 1
            if 0 <= idx < len(order):
                ref_name = order[idx]
            else:
                print(f"Error: invalid number: {user_input}")
                sys.exit(1)
        else:
            ref_name = user_input

        if ref_name not in sequences:
            print(f"Error: reference '{ref_name}' not found.")
            sys.exit(1)
        print(f"Reference selected: {ref_name}")

    #Check alignment
    aln_length = check_alignment(sequences)
    print(f"\nAlignment length: {aln_length:,} columns")

    ref_seq = sequences[ref_name]
    ref_non_gap = sum(1 for b in ref_seq if b != "-")
    print(f"Reference positions (excluding gaps): {ref_non_gap:,}\n")

    #Compute conservation
    results = []
    for name in order:
        if name == ref_name:
            continue
        pct, compared, identical = compute_identity(ref_seq, sequences[name])
        results.append({
            "species": name,
            "compared_positions": compared,
            "identical_positions": identical,
            "pct_identity": round(pct, 4),
        })
    mean = sum(r["pct_identity"] for r in results) / len(results) if results else 0.0

    #Display
    print(f"{'Species':<50} {'% Identity':>12} {'Identical':>12} {'Compared':>12}")
    for r in sorted(results, key=lambda x: x["pct_identity"], reverse=True):
        print(f"{r['species']:<50} {r['pct_identity']:>11.2f} {r['identical_positions']:>12,} {r['compared_positions']:>12,}")

    #Export CSV
    out_csv = Path("/Xnfs/lbmcdb/Semon_team/jcaron/id_scores/") / (fasta_path.stem + "_conservation.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["species", "pct_identity", "identical_positions", "compared_positions"])
        writer.writeheader()
        # Reference as first row
        writer.writerow({
            "species": ref_name + " (reference)",
            "pct_identity": 100.0,
            "identical_positions": ref_non_gap,
            "compared_positions": ref_non_gap,
        })
        for r in sorted(results, key=lambda x: x["pct_identity"], reverse=True):
            writer.writerow(r)
        writer.writerow({
            "species": "MEAN",
            "pct_identity": round(mean, 4),
            "identical_positions": "",
            "compared_positions": "",
        })
    print(f"\nCSV exported: {out_csv}")
    print(f"Mean % identity (vs reference): {mean:.2f}%\n")


if __name__ == "__main__":
    main()

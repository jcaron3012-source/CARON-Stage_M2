#!/usr/bin/env python3

import os
import re
import argparse
from collections import defaultdict

'''Merges FASTA fragments based on a BED sequence. Returns a FASTA
file named {gene}_{annotation}_{n_species_aligned}.fa'''

def parse_fasta(filepath):
    sequences = {}
    current_sp = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_sp = line[1:].split()[0]
                sequences[current_sp] = []
            elif current_sp:
                sequences[current_sp].append(line)
    return {sp: "".join(seqs) for sp, seqs in sequences.items()}


def parse_filename(filename):
    basename = os.path.basename(filename)
    m = re.search(r'_(chr[^_]+)_(\d+)_(\d+)\.fa$', basename)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return None


def sanitize(value):
    '''Cleanup of alignments'''
    value = value.strip()
    if not value or value in (".", "-"):
        return "unknown"
    value = re.sub(r'[^\w.-]+', '_', value)
    return value.strip('_') or "unknown"


def parse_bed_line(line):
    line = line.strip()
    if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
        return None
    parts = re.split(r'[\t ]+', line)
    if len(parts) < 8:
        return None
    try:
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        gene = sanitize(parts[3])
        annotation = sanitize(parts[7])
        if chrom.startswith("chr") and start >= 0 and end > start:
            return (chrom, start, end, gene, annotation)
    except ValueError:
        pass
    return None


def load_bed(bed_file, debug=False):
    regions = []
    failed = []
    with open(bed_file) as f:
        for i, line in enumerate(f):
            parsed = parse_bed_line(line)
            if parsed:
                regions.append(parsed)
                if debug and len(regions) <= 5:
                    print(f"line {i + 1}: {parsed}")
            else:
                if line.strip():
                    failed.append((i + 1, line.strip()[:80]))
    if failed:
        print(f"{len(failed)} BED lines not parsed")
        for lineno, content in failed[:3]:
            print(f"{lineno}: {content!r}")
    return regions


def get_unique_path(outdir, base_name, ext=".fa"):
    candidate = os.path.join(outdir, f"{base_name}{ext}")
    if not os.path.exists(candidate):
        return candidate
    i = 2
    while True:
        candidate = os.path.join(outdir, f"{base_name}_{i}{ext}")
        if not os.path.exists(candidate):
            return candidate
        i += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_dir", required=True)
    parser.add_argument("--maf_file", required=True)
    parser.add_argument("--bed_file", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    maf_file = args.maf_file
    bed_file = args.bed_file
    print(f"Processing {args.fasta_dir}, with {maf_file} and {bed_file}...")
    if not os.path.exists(bed_file):
        print(f"WARNING: No BED file found for {maf_file}, skipping.")
        return

    fa_by_chrom = defaultdict(list)

    if not os.path.exists(args.fasta_dir):
        print(f"WARNING: No FASTA directory found for {maf_file}, skipping.")
        return

    n_indexed = 0
    for fname in os.listdir(args.fasta_dir):
        coords = parse_filename(fname)
        if coords:
            path = os.path.join(args.fasta_dir, fname)
            chrom, s, e = coords
            fa_by_chrom[chrom].append((s, e, path))
            n_indexed += 1

    print(f"{n_indexed} .fa files indexed for {maf_file}.")

    bed_regions = load_bed(bed_file, debug=args.debug)
    print(f"{len(bed_regions)} BED regions loaded for {maf_file}.")

    if not bed_regions:
        print(f"FATAL: No BED regions found for {maf_file}, launch with --debug.")
        return

    merged_count = 0
    skipped_count = 0

    for (bed_chr, bed_start, bed_end, gene, annotation) in bed_regions:
        candidates = fa_by_chrom.get(bed_chr, [])
        blocks = sorted(
            [(s, e, path) for (s, e, path) in candidates if s >= bed_start and e <= bed_end],
            key=lambda x: x[0]
        )

        if not blocks:
            print(f"FATAL: No block found for {bed_chr}_{bed_start}_{bed_end} in {maf_file}")
            skipped_count += 1
            continue

        all_species = set()
        block_data = []
        for (blk_s, blk_e, path) in blocks:
            fasta = parse_fasta(path)
            block_data.append((blk_s, blk_e, fasta))
            all_species.update(fasta.keys())

        merged = defaultdict(str)
        for (blk_s, blk_e, fasta) in block_data:
            blk_aln_len = len(next(iter(fasta.values()))) if fasta else 0
            for sp in all_species:
                if sp in fasta:
                    merged[sp] += fasta[sp]
                else:
                    merged[sp] += "-" * blk_aln_len

        n_sp = len(all_species)

        base_name = f"{gene}_{annotation}_{n_sp}"
        out_path = get_unique_path(args.outdir, base_name)
        out_fname = os.path.basename(out_path)

        with open(out_path, "w") as out:
            for sp in sorted(merged.keys()):
                out.write(f">{sp}\n")
                seq = merged[sp]
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")
        merged_count += 1
        if merged_count % 100 == 0:
            print(f"{merged_count} treated for {maf_file}")
        print(f"Treated {out_fname}")

    print(f"\nDone for {maf_file}. {merged_count} files, {skipped_count} regions without blocks.")


if __name__ == "__main__":
    main()

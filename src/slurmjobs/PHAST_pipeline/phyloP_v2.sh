#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=Lake
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --job-name=bedtoolsmap
#SBATCH --output=/Xnfs/lbmcdb/Semon_team/jcaron/log/bedtoolsmap.out
#SBATCH --error=/Xnfs/lbmcdb/Semon_team/jcaron/log/bedtoolsmap.err
set -euo pipefail

# Calculates the mean phyloP score for each CRE

START=$(date +%s)

MUS="/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/mus_links.bed"
PHAST_DIR="/Xnfs/lbmcdb/Semon_team/jcaron/files/mdongScores/PhyloP_53"
OUT="/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/mus_phyloP_v2.bed"

# Scratch directory: use the SLURM-provided local tmpdir if available,
# otherwise fall back to a per-user/per-job directory under /tmp.
# Always cleaned up on exit.
LOCAL_TMP=${SLURM_TMPDIR:-/tmp/${USER}_$SLURM_JOB_ID}
mkdir -p "$LOCAL_TMP"
trap 'rm -rf "$LOCAL_TMP"' EXIT

echo "[1/4] Sort of $MUS..."
sort -k1,1 -k2,2n --buffer-size=2G --temporary-directory="$LOCAL_TMP" "$MUS" > "$LOCAL_TMP/mus_sorted.bed"

# List of chromosomes present in the CRE file; each one is processed independently
CHROMS=$(cut -f1 "$LOCAL_TMP/mus_sorted.bed" | sort -u)
echo "Chromosomes: $(echo "$CHROMS" | wc -l)"

echo "[2/4] Indexing phyloP window files..."
# Build a lookup table of the per-window phyloP files: chrom, window start,
# window end, and the file path. File names encode the chromosome, window
# start, and window size, e.g. "<chrom>.<start>.<size>_scoresPhyloP_250.wig.bed".
# Index: chrom <tab> winstart <tab> winend <tab> filepath
> "$LOCAL_TMP/phast_index.tsv"
for f in "$PHAST_DIR"/*_scoresPhyloP_250.wig.bed; do
    base=$(basename "$f")
    prefix="${base%_scoresPhyloP_250.wig.bed}"
    chrom="${prefix%%.*}"
    rest="${prefix#*.}"
    wstart="${rest%%.*}"
    size="${rest#*.}"
    wend=$(( wstart + size ))
    printf "%s\t%s\t%s\t%s\n" "$chrom" "$wstart" "$wend" "$f"
done > "$LOCAL_TMP/phast_index.tsv"
# Sort the index by chrom + window start so files are streamed in genomic order
sort -k1,1 -k2,2n "$LOCAL_TMP/phast_index.tsv" -o "$LOCAL_TMP/phast_index.tsv"
echo "Indexed $(wc -l < "$LOCAL_TMP/phast_index.tsv") PhyloP window files."

# awk script: accumulates mean PhyloP score per feature via streaming,
# without ever sorting or fully loading the PhyloP window files.
#
# Approach: features are read first (NR == FNR) and indexed both by start
# and by end position. As phyloP score lines are then streamed in genomic
# order, a running sum/count of scores is maintained; each feature snapshots
# that running total the moment the stream reaches its start, and again the
# moment it reaches its end. The mean for a feature is simply the difference
# of those two snapshots divided by the difference in counts -- this avoids
# ever storing per-base scores or re-scanning the phyloP stream per feature.
cat > "$LOCAL_TMP/map_mean.awk" << 'AWKEOF'
BEGIN {
    FS = "\t"; OFS = "\t"
    runsum = 0
    runcnt = 0
    sorted = 0
}

NR == FNR {
    n++
    fstart[n] = $2 + 0
    fend[n]   = $3 + 0
    line[n]   = $0
    next
}

# Build sort-by-start and sort-by-end index arrays once, right before
# processing the first streamed phyloP line. Insertion sort is fine:
# n is at most a few thousand features per chromosome.
!sorted {
    for (k = 1; k <= n; k++) { byStart[k] = k; byEnd[k] = k }
    for (a = 2; a <= n; a++) {
        keyS = byStart[a]; vS = fstart[keyS]
        b = a - 1
        while (b >= 1 && fstart[byStart[b]] > vS) { byStart[b+1] = byStart[b]; b-- }
        byStart[b+1] = keyS

        keyE = byEnd[a]; vE = fend[keyE]
        b = a - 1
        while (b >= 1 && fend[byEnd[b]] > vE) { byEnd[b+1] = byEnd[b]; b-- }
        byEnd[b+1] = keyE
    }
    startPtr = 1
    endPtr = 1
    sorted = 1
}

{
    pos = $2 + 0
    score = $5 + 0

    # Snapshot running totals for every feature whose start we've now reached.
    while (startPtr <= n && fstart[byStart[startPtr]] <= pos) {
        k = byStart[startPtr]
        sum_at_start[k] = runsum
        cnt_at_start[k] = runcnt
        startPtr++
    }

    # Finalize every feature whose end we've now reached.
    while (endPtr <= n && fend[byEnd[endPtr]] <= pos) {
        k = byEnd[endPtr]
        sum_at_end[k] = runsum
        cnt_at_end[k] = runcnt
        endPtr++
    }

    runsum += score
    runcnt += 1
}

END {
    # Edge case: no phyloP lines were streamed at all for this chromosome.
    if (!sorted) {
        for (k = 1; k <= n; k++) { byEnd[k] = k }
        for (a = 2; a <= n; a++) {
            keyE = byEnd[a]; vE = fend[keyE]
            b = a - 1
            while (b >= 1 && fend[byEnd[b]] > vE) { byEnd[b+1] = byEnd[b]; b-- }
            byEnd[b+1] = keyE
        }
        endPtr = 1
    }
    # Finalize any feature whose end lies beyond the last base streamed.
    while (endPtr <= n) {
        k = byEnd[endPtr]
        sum_at_end[k] = runsum
        cnt_at_end[k] = runcnt
        endPtr++
    }
    for (i = 1; i <= n; i++) {
        c = cnt_at_end[i] - cnt_at_start[i]
        s = sum_at_end[i] - sum_at_start[i]
        if (c > 0) print line[i], s / c; else print line[i], "."
    }
}
AWKEOF

# Processes a single chromosome: extracts its CREs and phyloP window files,
# then streams the window files through the awk script above to compute
# each CRE's mean phyloP score. Run in parallel across chromosomes below.
run_chrom() {
    local CHR=$1
    local TMPDIR=$2
    local INDEX=$3
    echo "--- Treating $CHR ---"

    awk -v chr="$CHR" '$1 == chr' "$TMPDIR/mus_sorted.bed" > "$TMPDIR/mus_${CHR}.bed"

    if [[ ! -s "$TMPDIR/mus_${CHR}.bed" ]]; then
        return
    fi

    awk -v chr="$CHR" '$1 == chr {print $4}' "$INDEX" > "$TMPDIR/files_${CHR}.txt"

    if [[ ! -s "$TMPDIR/files_${CHR}.txt" ]]; then
        echo "WARN: No PhyloP window files for $CHR" >&2
        awk 'BEGIN{OFS="\t"} {print $0, "."}' "$TMPDIR/mus_${CHR}.bed" > "$TMPDIR/out_phyloP_${CHR}.bed"
        rm -f "$TMPDIR/mus_${CHR}.bed" "$TMPDIR/files_${CHR}.txt"
        return
    fi

    echo "Streaming $(wc -l < "$TMPDIR/files_${CHR}.txt") PhyloP window file(s) for $CHR..."
    # Stream all window files for this chromosome through awk, never
    # materializing them concatenated or sorted on disk/RAM.
    xargs -a "$TMPDIR/files_${CHR}.txt" cat \
        | awk -f "$TMPDIR/map_mean.awk" "$TMPDIR/mus_${CHR}.bed" - \
        > "$TMPDIR/out_phyloP_${CHR}.bed"

    rm -f "$TMPDIR/mus_${CHR}.bed" "$TMPDIR/files_${CHR}.txt"
}
export -f run_chrom

echo "[3/4] Streaming PhyloP mean per chromosome..."
# Fan out across chromosomes with GNU parallel (4 concurrent jobs); abort
# the whole run as soon as any chromosome fails.
echo "$CHROMS" | parallel --halt now,fail=1 -j 4 run_chrom {} "$LOCAL_TMP" "$LOCAL_TMP/phast_index.tsv"

echo "[4/4] Final assembly..."
# Concatenate the per-chromosome results back together and re-sort into a
# single genome-wide, coordinate-sorted output file
find "$LOCAL_TMP" -maxdepth 1 -name "out_phyloP_*.bed" -print0 \
    | xargs -0 cat \
    | sort -k1,1 -k2,2n --buffer-size=2G --temporary-directory="$LOCAL_TMP" \
    > "$LOCAL_TMP/result.bed"

echo "$(wc -l < "$LOCAL_TMP/result.bed") lines in $OUT"
cp "$LOCAL_TMP/result.bed" "$OUT"

END=$(date +%s); ELAPSED=$(( END - START ))
printf "Total time: %02d:%02d:%02d:%02d (jj:hh:mm:ss)\n" \
  $(( ELAPSED/86400 )) $(( (ELAPSED%86400)/3600 )) \
  $(( (ELAPSED%3600)/60 )) $(( ELAPSED%60 ))

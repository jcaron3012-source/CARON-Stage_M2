#!/usr/bin/env bash
set -uo pipefail

#Takes as input the results of hal2maf (MAF chunks), converts every alignment in FASTA format and concatenates fragmented sequences into one single cis-regulatory regions with the information provided on the BED file.

DIR="/Xnfs/lbmcdb/Semon_team/jcaron"
MAFDIR="$DIR/files/maffiles"
BEDDIR="$DIR/files/bedfiles/chunks"
OUTDIR="$DIR/files/fasta_msa"
MERGE_SCRIPT="$DIR/src/python_scripts/merge_msa_by_region.py"

TMPDIR_BLOCKS="$TMPDIR/blocks"
MARKERDIR="$OUTDIR/.markers"       
STAGING_PARENT="$OUTDIR/.pending" 

echo "Writing blocks in $TMPDIR_BLOCKS ..."
mkdir -p "$TMPDIR_BLOCKS" || { echo "FATAL: can't create $TMPDIR_BLOCKS"; exit 1; }
mkdir -p "$OUTDIR" || { echo "FATAL: can't create $OUTDIR"; exit 1; }
mkdir -p "$MARKERDIR" || { echo "FATAL: can't create $MARKERDIR"; exit 1; }
mkdir -p "$STAGING_PARENT" || { echo "FATAL: can't create $STAGING_PARENT"; exit 1; }

move_chunk_outputs() {
    local staging="$1"
    local dest="$2"
    local f base candidate i
    for f in "$staging"/*.fa; do
        [[ -e "$f" ]] || continue
        base=$(basename "$f" .fa)
        candidate="$dest/${base}.fa"
        i=2
        while [[ -e "$candidate" ]]; do
            candidate="$dest/${base}_${i}.fa"
            (( i++ ))
        done
        mv "$f" "$candidate"
    done
}

# ==============================================================================
write_block() {
    #Writes the FASTA files for each block
    local block="$1"
    local blockdir="$2"
    local n_species
    n_species=$(grep -c '^s' <<< "$block")
    local sline
    sline=$(grep '^s' <<< "$block" | head -1 | tr '\t' ' ')
    if [[ -z "$sline" ]]; then
        echo "WARNING: no 's' line found in block, skipping" >&2
        return 1
    fi
    local _ src start size
    read -r _ src start size _ _ _ <<< "$sline"
    if [[ -z "$src" ]]; then
        echo "WARNING: could not parse 's' line: $sline" >&2
        return 1
    fi
    local species="${src%%.*}"
    local chrom="${src#*.}"
    local end=$(( start + size ))
    local outfile="${blockdir}/${n_species}-${species}_${chrom}_${start}_${end}.fa"
    printf '##maf version=1\n\n%s\n' "$block" > /tmp/_block.maf
    if ! msa_view /tmp/_block.maf --in-format MAF --out-format FASTA > "$outfile" 9>&-; then
        echo "WARNING: msa_view failed on block -> $outfile" >&2
        rm -f "$outfile"
        return 1
    fi
}

# ==============================================================================
print_progress() {
    local done=$1
    local total=$2
    local pct=$(( total > 0 ? done * 100 / total : 100 ))
    local filled=$(( pct / 2 ))
    local empty=$(( 50 - filled ))
    local bar="["
    bar+=$(printf '#%.0s' $(seq 1 $filled))
    bar+=$(printf ' %.0s' $(seq 1 $empty))
    bar+="]"
    echo "$bar $pct% ($done/$total)"
}

for MAF in "$MAFDIR"/chunk_*.maf; do
    [[ -e "$MAF" ]] || continue
    MAF_BASENAME=$(basename "$MAF" .maf)
    BED_FILE="$BEDDIR/${MAF_BASENAME}.bed"
    BLOCKDIR="$TMPDIR_BLOCKS/$MAF_BASENAME"
    STAGING="$STAGING_PARENT/$MAF_BASENAME"
    MARKER="$MARKERDIR/${MAF_BASENAME}.done"
    
    #Check if the file is already done
    if [[ -f "$MARKER" ]]; then
        echo "$MAF_BASENAME already done"
        continue
    fi

    if [[ ! -f "$BED_FILE" ]]; then
        echo "FATAL: No BED file found for $MAF_BASENAME (expected: $BED_FILE)"
        continue
    fi

    rm -rf "$BLOCKDIR" "$STAGING"
    mkdir -p "$BLOCKDIR" "$STAGING"

    echo ""
    echo "=== $MAF_BASENAME ==="
    echo "Step 1: MAF to FASTA"

    total_alignments=$(grep -c '^a' "$MAF")
    echo "Expected alignments: $total_alignments"

    files_created=0
    files_skipped=0
    tmp=""
    
    #Main command: FASTA conversion
    while IFS= read -r -u9 line || [[ -n "$line" ]]; do
        if [[ "$line" =~ ^a ]]; then
            if [[ -n "$tmp" ]]; then
                if ( write_block "$tmp" "$BLOCKDIR" ); then
                    (( files_created++ ))
                else
                    (( files_skipped++ ))
                fi
                print_progress $(( files_created + files_skipped )) "$total_alignments"
            fi
            tmp="$line"
        elif [[ "$line" =~ ^[sieq] ]]; then
            tmp+=$'\n'"$line"
        fi
    done 9< "$MAF"

    #Last file
    if [[ -n "$tmp" ]]; then
        if write_block "$tmp" "$BLOCKDIR"; then
            (( files_created++ ))
        else
            (( files_skipped++ ))
        fi
        print_progress $(( files_created + files_skipped )) "$total_alignments"
    fi
    echo "$files_created blocks created for $MAF_BASENAME in $BLOCKDIR ($files_skipped ignored)"

    echo ""
    echo "Step 2: fusion by BED region"
    echo "Processing $MAF_BASENAME in $BLOCKDIR with BED: $BED_FILE"

    #Merging the fragmented outputs
    python3 "$MERGE_SCRIPT" \
        --fasta_dir "$BLOCKDIR" \
        --maf_file "$MAF" \
        --bed_file "$BED_FILE" \
        --outdir "$STAGING"

    if [[ $? -ne 0 ]]; then
        echo "FATAL: python script failed for $MAF_BASENAME."
        echo "Blocks conserved in $BLOCKDIR and $STAGING for inspection."
        exit 1
    fi

    echo "Step 3: Checking outputs $MAF_BASENAME in $OUTDIR"
    move_chunk_outputs "$STAGING" "$OUTDIR"
    rm -rf "$STAGING" "$BLOCKDIR"

    touch "$MARKER"
done

rm -rf "$TMPDIR_BLOCKS" "$STAGING_PARENT"

# ==============================================================================
echo ""
echo "Step 4: Removing alignments with only two species"

removed=0
for f in "$OUTDIR"/*.fa; do
    [[ -e "$f" ]] || continue
    n=$(grep -c '^>' "$f")
    if [[ "$n" -le 2 ]]; then
        rm "$f"
        (( removed++ ))
    fi
done
echo "Done: $removed files removed. $(ls "$OUTDIR" | wc -l) files remaining in $OUTDIR"

# ==============================================================================
echo ""
echo "Step 5: replacing N by -"
for f in "$OUTDIR"/*.fa; do
    sed -i 's/N/-/g' "$f"
    sed -i 's/-annospalax_galili/Nannospalax_galili/g' "$f"
done
echo "msa_view done."

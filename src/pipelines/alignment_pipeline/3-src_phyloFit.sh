#!/usr/bin/env bash

#Adjusts a REV model on all FASTA files, needs a Newick format file for tree topology. 
#Returns a tree for each sequence.

OUTDIR="/Xnfs/lbmcdb/Semon_team/jcaron/files/modfiles"
TREE="/Xnfs/lbmcdb/Semon_team/jcaron/files/nhfiles/rodents_noamericans.nh"
mkdir -p "$OUTDIR"

fa=$1

if [[ ! -f "$fa" ]]; then
    echo "FATAL: $fa does not exist." >&2
    exit 1
fi

base=$(basename "${fa%.fa}")
mod="${OUTDIR}/${base}.mod"

echo "Processing: $base"

if grep -q "^TREE:" "$mod" 2>/dev/null; then
    echo "Already done, skipping."
    exit 0
fi

#Main command
phyloFit --tree "$TREE" --subst-mod REV --out-root "${OUTDIR}/${base}" "$fa" || true

#Retry if it failed
if ! grep -q "^TREE:" "$mod" 2>/dev/null; then
    echo "First attempt failed, retrying..."
    phyloFit --tree "$TREE" --subst-mod REV --out-root "${OUTDIR}/${base}" "$fa" || true
fi

if grep -q "^TREE:" "$mod" 2>/dev/null; then
    echo "Success."
else
    echo "FAILED: $base" >&2
    exit 1
fi

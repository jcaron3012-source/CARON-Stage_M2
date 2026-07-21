#!/bin/bash
#SBATCH --job-name=phyloFit_summary
#SBATCH --output=/Xnfs/lbmcdb/Semon_team/jcaron/log/phyloFit_summary.out
#SBATCH --error=/Xnfs/lbmcdb/Semon_team/jcaron/log/phyloFit_summary.err
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --partition=Lake

#Do not run, use 3_phyloFit instead !!

LOG_DIR=/Xnfs/lbmcdb/Semon_team/jcaron/log
FASTA_DIR=/Xnfs/lbmcdb/Semon_team/jcaron/files/fasta_msa/
OUTPUT_DIR=/Xnfs/lbmcdb/Semon_team/jcaron/files/modfiles/
IS_RERUN=${1:-0}

# Compute and display total elapsed time
START=$(cat "${LOG_DIR}/phyloFit_START.tmp")
END=$(date +%s)
ELAPSED=$(( END - START ))
DAYS=$(( ELAPSED / 86400 ))
HOURS=$(( (ELAPSED % 86400) / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))
printf "Total time : %02d:%02d:%02d:%02d (dd:hh:mm:ss)\n" $DAYS $HOURS $MINS $SECS

# Count input and output files
N_INPUT=$(ls "$FASTA_DIR" | wc -l)
N_OUTPUT=$(ls "$OUTPUT_DIR" | wc -l)
echo "Input files  : $N_INPUT"
echo "Output files : $N_OUTPUT"

if [ "$N_OUTPUT" -lt "$N_INPUT" ]; then
    MISSING=$(( N_INPUT - N_OUTPUT ))
    echo "WARNING: $MISSING missing output file(s)."
    if [ "$IS_RERUN" -eq 0 ]; then
        echo "Triggering a single automatic rerun for missing files..."
        # Identify missing indices by comparing input filenames to output files
        MISSING_INDICES=()
        i=1
        for f in $(ls "$FASTA_DIR"); do
            BASENAME=$(basename "$f" .fa)
            if [ ! -f "${OUTPUT_DIR}${BASENAME}.mod" ]; then
                MISSING_INDICES+=($i)
            fi
            i=$(( i + 1 ))
        done
        ARRAY_STR=$(IFS=,; echo "${MISSING_INDICES[*]}")
        echo "Re-submitting indices : $ARRAY_STR"

        # Record new start time for the rerun
        date +%s > "${LOG_DIR}/phyloFit_START.tmp"

        JOB_ID=$(sbatch --parsable --array=${ARRAY_STR}%50 /Xnfs/lbmcdb/Semon_team/jcaron/src/slurmjobs/alignment_pipeline/3_phyloFit_command.sh 1)
        echo "Rerun job submitted : $JOB_ID"

        # Submit summary again after rerun, with IS_RERUN=1 to prevent further reruns
        sbatch --dependency=afterok:"$JOB_ID" /Xnfs/lbmcdb/Semon_team/jcaron/src/slurmjobs/alignment_pipeline/3_phyloFit_summary.sh 1
    else
        echo "ERROR: Rerun already attempted — files still missing. Manual inspection required."
    fi
else
    echo "SUCCESS: All output files are present."
    rm -f "${LOG_DIR}/phyloFit_START.tmp"
fi

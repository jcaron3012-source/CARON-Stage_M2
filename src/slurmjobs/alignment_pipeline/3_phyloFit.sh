#!/bin/bash

#Submit this job for phyloFit !

FASTA_DIR="/Xnfs/lbmcdb/Semon_team/jcaron/files/fasta_msa"
LOG_DIR="/Xnfs/lbmcdb/Semon_team/jcaron/log"
LIST_FILE="${LOG_DIR}/file_list.txt"
mkdir -p "$LOG_DIR"

# File list to check if we did them all
find "$FASTA_DIR" -maxdepth 1 -type f > "$LIST_FILE"

# Array of the number of files
N=$(wc -l < "$LIST_FILE")

MAX_ARRAY_SIZE=$(scontrol show config 2>/dev/null | awk -F= '/MaxArraySize/{gsub(/ /,"",$2); print $2}')
MAX_ARRAY_SIZE=${MAX_ARRAY_SIZE:-1000}

FILES_PER_JOB=$(( (N + MAX_ARRAY_SIZE - 2) / (MAX_ARRAY_SIZE - 1) ))
[[ "$FILES_PER_JOB" -lt 1 ]] && FILES_PER_JOB=1

NUM_TASKS=$(( (N + FILES_PER_JOB - 1) / FILES_PER_JOB ))
echo "Number of files : $N"
echo "MaxArraySize : $MAX_ARRAY_SIZE"
echo "Slurm jobs (groups of $FILES_PER_JOB) : $NUM_TASKS"

# Start time
date +%s > "${LOG_DIR}/phyloFit_START.tmp"

# Job array submit
JOB_ID=$(sbatch --parsable --array=1-${NUM_TASKS}%50 /Xnfs/lbmcdb/Semon_team/jcaron/src/slurmjobs/alignment_pipeline/3_phyloFit_command.sh "$FILES_PER_JOB")

# Checking if the submission is done right
if [ -n "$JOB_ID" ] && [[ "$JOB_ID" =~ ^[0-9]+ ]]; then
    echo "Job array submitted with success : $JOB_ID"

    # Summary job
    sbatch --dependency=afterok:"$JOB_ID" /Xnfs/lbmcdb/Semon_team/jcaron/src/slurmjobs/alignment_pipeline/3_phyloFit_summary.sh 0
    echo "Submitting summary job, waiting for array to finish."
else
    echo "[ERROR] Job array failed."
    exit 1
fi

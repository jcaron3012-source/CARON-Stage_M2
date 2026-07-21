#!/bin/bash
#SBATCH --job-name=phyloFit_job
#SBATCH --output=/Xnfs/lbmcdb/Semon_team/jcaron/log/phyloFit_%a.out
#SBATCH --error=/Xnfs/lbmcdb/Semon_team/jcaron/log/phyloFit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=Cascade

#Do not run, use 3_phyloFit.sh !!

#Activating PHAST environment
source /Xnfs/lbmcdb/Semon_team/jcaron/utils/miniconda3/etc/profile.d/conda.sh
conda activate /Xnfs/lbmcdb/Semon_team/jcaron/utils/phast_env

BATCH_SIZE=$1
LIST_FILE="/Xnfs/lbmcdb/Semon_team/jcaron/log/file_list.txt"

# Number of lines to do
START_LINE=$(( ($SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END_LINE=$(( $SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
echo "Starting job ${SLURM_ARRAY_TASK_ID}"
echo "From line $START_LINE to $END_LINE of $LIST_FILE"

sed -n "${START_LINE},${END_LINE}p" "$LIST_FILE" | while read -r FILE_PATH; do
    if [ -n "$FILE_PATH" ] && [ -f "$FILE_PATH" ]; then
        echo "Processing: $FILE_PATH"
        START_TIME=$(date +%s)
        # starting the phyloFit pipeline
        bash /Xnfs/lbmcdb/Semon_team/jcaron/src/pipelines/alignment_pipeline/3-src_phyloFit.sh "$FILE_PATH"
        END_TIME=$(date +%s)
        ELAPSED=$(( END_TIME - START_TIME ))
        echo "File done in $ELAPSED seconds."
    else
        echo "[INFO] $FILE_PATH not found, onto the next."
    fi
done
echo "End of command."

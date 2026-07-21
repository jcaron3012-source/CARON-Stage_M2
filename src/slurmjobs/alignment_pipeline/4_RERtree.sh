#!/bin/bash
#SBATCH --job-name=RERtree_job
#SBATCH --output=log/RERtree.out
#SBATCH --error=log/RERtree.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

START=$(date +%s)

bash /Xnfs/lbmcdb/Semon_team/jcaron/pipelines/RERtree.sh

END=$(date +%s)
ELAPSED=$(( END - START ))

DAYS=$(( ELAPSED / 86400 ))
HOURS=$(( (ELAPSED % 86400) / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

printf "Time taken : %02d:%02d:%02d:%02d (jj:hh:mm:ss)\n" $DAYS $HOURS $MINS $SECS


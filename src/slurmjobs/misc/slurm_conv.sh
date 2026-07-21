#!/bin/bash
#SBATCH --job-name=conv_job
#SBATCH --output=log/conv.out
#SBATCH --error=log/conv.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

START=$(date +%s)

bash /Xnfs/lbmcdb/Semon_team/jcaron/src/pipelines/misc/conversion.sh

END=$(date +%s)
ELAPSED=$(( END - START ))

DAYS=$(( ELAPSED / 86400 ))
HOURS=$(( (ELAPSED % 86400) / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

printf "Time: %02d:%02d:%02d:%02d (jj:hh:mm:ss)\n" $DAYS $HOURS $MINS $SECS

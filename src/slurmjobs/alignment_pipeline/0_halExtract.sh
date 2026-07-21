#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --partition=Lake
#SBATCH --mem=150G
#SBATCH --time=08:00:00
#SBATCH --job-name=halExtract
#SBATCH --output=/Xnfs/lbmcdb/Semon_team/jcaron/log/halExtract.out
#SBATCH --error=/Xnfs/lbmcdb/Semon_team/jcaron/log/halExtract.err

START=$(date +%s)

/Xnfs/abc/charliecloud_bin/ch-run --write --bind=/scratch/:/scratch/ --bind=/Xnfs:/Xnfs /Xnfs/abc/charliecloud/img/cactus_v3.1.4.sqfs -- bash /Xnfs/lbmcdb/Semon_team/jcaron/src/pipelines/alignment_pipeline/0_src_halExtract.sh

END=$(date +%s)
ELAPSED=$(( END - START ))

DAYS=$(( ELAPSED / 86400 ))
HOURS=$(( (ELAPSED % 86400) / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

printf "Time taken : %02d:%02d:%02d:%02d (jj:hh:mm:ss)\n" $DAYS $HOURS $MINS $SECS


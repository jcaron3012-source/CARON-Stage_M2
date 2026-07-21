#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=Lake
#SBATCH --mem-per-cpu=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=hal2maf
#SBATCH --output=/Xnfs/lbmcdb/Semon_team/jcaron/log/hal2fasta.out
#SBATCH --error=/Xnfs/lbmcdb/Semon_team/jcaron/log/hal2fasta.err

START=$(date +%s)

/Xnfs/abc/charliecloud_bin/ch-run --write --bind=/scratch/:/scratch/ --bind=/Xnfs:/Xnfs /Xnfs/abc/charliecloud/img/cactus_v3.1.4.sqfs -- bash /Xnfs/lbmcdb/Semon_team/jcaron/src/pipelines/misc/hal2fasta.sh

END=$(date +%s)
ELAPSED=$(( END - START ))

DAYS=$(( ELAPSED / 86400 ))
HOURS=$(( (ELAPSED % 86400) / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

printf "Time taken : %02d:%02d:%02d:%02d (jj:hh:mm:ss)\n" $DAYS $HOURS $MINS $SECS

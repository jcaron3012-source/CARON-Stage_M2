#!/bin/bash
#SBATCH --job-name=msa_view_job
#SBATCH --output=log/msa_view.out
#SBATCH --error=log/msa_view.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p Lake
#SBATCH --time=01-00:00:00
#SBATCH --mem=4G

#Activating PHAST environment
source /Xnfs/lbmcdb/Semon_team/jcaron/utils/miniconda3/bin/activate
conda activate /Xnfs/lbmcdb/Semon_team/jcaron/utils/phast_env

START=$(date +%s)

bash /Xnfs/lbmcdb/Semon_team/jcaron/src/pipelines/alignment_pipeline/2-src_msa_view.sh

END=$(date +%s)
ELAPSED=$(( END - START ))

DAYS=$(( ELAPSED / 86400 ))
HOURS=$(( (ELAPSED % 86400) / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

printf "Time : %02d:%02d:%02d:%02d (jj:hh:mm:ss)\n" $DAYS $HOURS $MINS $SECS

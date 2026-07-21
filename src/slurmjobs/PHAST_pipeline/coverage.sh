#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=Lake
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --job-name=coverage
#SBATCH --output=/Xnfs/lbmcdb/Semon_team/jcaron/log/coverage.out
#SBATCH --error=/Xnfs/lbmcdb/Semon_team/jcaron/log/coverage.err
set -euo pipefail

#Intersects the regions predicted by Viterbi and a BED file

START=$(date +%s)

CRE="/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/mus_phastcons_v2.bed"
VTB="/Xnfs/lbmcdb/Semon_team/jcaron/files/mdongScores/PhastConsRegions_mouse_onlychr_v3_mdong_10Mb_53RODENTS_MTDF.bed"
OUT="/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/coverage_CREs.bed"

bedtools coverage -a $CRE -b $VTB > $OUT

echo "Viterbi: "
head -5 $VTB
echo "-----------------------------------------------------------------------------------"
echo "Out: "
head -5 $OUT

END=$(date +%s); ELAPSED=$(( END - START ))
printf "Total time: %02d:%02d:%02d:%02d (jj:hh:mm:ss)\n" \
  $(( ELAPSED/86400 )) $(( (ELAPSED%86400)/3600 )) \
  $(( (ELAPSED%3600)/60 )) $(( ELAPSED%60 ))

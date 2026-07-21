#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=halLiftOver
#SBATCH --output=/Xnfs/abc/nf_scratch/mrouyer/halper/log/halLiftover.out
#SBATCH --error=/Xnfs/abc/nf_scratch/mrouyer/halper/log/halLiftover.err
#SBATCH --time=04:00:00

halLiftover \
/Xnfs/lbmcdb/Semon_team/jcaron/files/halfiles/rodents_noamericans.hal \
Mesocricetus_auratus \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1.bed \
Mus_musculus \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38.bed \
--bedType 3

halLiftover \
/Xnfs/lbmcdb/Semon_team/jcaron/files/halfiles/rodents_noamericans.hal \
Mesocricetus_auratus \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_summits_no_colnames_lift_nflo_chainfile_to_maur1.bed \
Mus_musculus \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_summits_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38.bed \
--bedType 3

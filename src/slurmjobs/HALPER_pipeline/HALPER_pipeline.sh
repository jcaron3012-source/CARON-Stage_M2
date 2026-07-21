#!/bin/bash

# Converts BED files of sequences and summits, and performs the HALPER
# pipeline to find orthologs between mouse (Mus musculus) and hamster
# (Mesocricetus auratus) regulatory regions.

##############################################################################
# 1. Add a peak_id column (chr-start-end) to each species' peak file
##############################################################################

python3 -c "  
with open('/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/HALPER_files/mus_by_molar_with_best_summit_median.bed') as f:  
    print(next(f).strip() + '\tpeak_id')  
    for line in f:  
        l = line.strip()  
        cols = l.split('\t')  
        print(f'{l}\t{cols[0]}-{cols[1]}-{cols[2]}')  
" > /Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id.bed

python3 -c "  
with open('/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/HALPER_files/ham_by_molar_with_best_summit_median.bed') as f:  
    print(next(f).strip() + '\tpeak_id')  
    for line in f:  
        l = line.strip()  
        cols = l.split('\t')  
        print(f'{l}\t{cols[0]}-{cols[1]}-{cols[2]}')  
" > /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id.bed

##############################################################################
# 2. Convert hamster chromosome names from RefSeq to GenBank accessions
##############################################################################

awk -F "\t" '
    
    # Mapping refseq to genebank
    NR==FNR {

        split($1, bits, " ")
        refseq = bits[1]
        genbank = bits[2]
        
        gsub(/\r/, "", refseq)
        gsub(/\r/, "", genbank)
        
        if (refseq != "") {
            map[refseq] = genbank
        }
        next
    }

    FNR==1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "summit_chr") {
                idx_summit = i
            }
        }
        print $0; next
    }


    {
        gsub(/\r/, "")
        
        if ($1 in map) {
            $1 = map[$1]
        }
        
        if (idx_summit && ($idx_summit in map)) {
            $idx_summit = map[$idx_summit]
        }
        
        print $0
    }
' OFS="\t" \
/Xnfs/lbmcdb/Semon_team/mrouyer/genomes/BCMMaur2/refseq_to_genbank.txt \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id.bed \
> /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names.bed


##############################################################################
# 3. Split each file into a peaks-only file and a summits-only file
##############################################################################

python3 -c "

file_input = '/Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id.bed'
peak_file_out = '/Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_peaks.bed'
summit_file_out = '/Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_summits.bed'

with open(file_input) as f, \
    open(peak_file_out, 'w') as f_peak, \
    open(summit_file_out, 'w') as f_summit:
    
    # header
    header = next(f).strip().split('\t')
    f_peak.write('\t'.join([header[0], header[1], header[2], header[3], header[4], header[-1]]) + '\n')
    f_summit.write('\t'.join(header[-7:]) + '\n')
    
    # data
    for line in f:
        cols = line.strip().split('\t')
        f_peak.write('\t'.join([cols[0], cols[1], cols[2], cols[3], cols[4], cols[-1]]) + '\n')
        f_summit.write('\t'.join(cols[-7:]) + '\n')
"

python3 -c "

file_input = '/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id.bed'
peak_file_out = '/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks.bed'
summit_file_out = '/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_summits.bed'

with open(file_input) as f, \
    open(peak_file_out, 'w') as f_peak, \
    open(summit_file_out, 'w') as f_summit:
    
    # header
    header = next(f).strip().split('\t')
    f_peak.write('\t'.join([header[0], header[1], header[2], header[3], header[4], header[-1]]) + '\n')
    f_summit.write('\t'.join(header[-7:]) + '\n')
    
    # data
    for line in f:
        cols = line.strip().split('\t')
        f_peak.write('\t'.join([cols[0], cols[1], cols[2], cols[3], cols[4], cols[-1]]) + '\n')
        f_summit.write('\t'.join(cols[-7:]) + '\n')
"

##############################################################################
# 4. Strip header rows from each peaks/summits file
##############################################################################

tail -n +2 /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_peaks.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_peaks_no_colnames.bed
tail -n +2 /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_summits.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_summits_no_colnames.bed

tail -n +2 /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames.bed
tail -n +2 /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_summits.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_summits_no_colnames.bed

tail -n +2 /Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_peaks.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_peaks_no_colnames.bed
tail -n +2 /Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_summits.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_summits_no_colnames.bed

# Strip the "chr" prefix from mouse chromosome names (peaks and summits)

sed -E 's/^chr([0-9]+)\t/\1\t/' mus_by_molar_with_best_summit_median_with_id_peaks_no_colnames.bed > mus_by_molar_with_best_summit_median_with_id_peaks_no_colnames_no_chr.bed
sed -E 's/^chr([0-9]+)\t/\1\t/' mus_by_molar_with_best_summit_median_with_id_summits_no_colnames.bed > mus_by_molar_with_best_summit_median_with_id_summits_no_colnames_no_chr.bed

##############################################################################
# 5. LiftOver: hamster BCMMaur2.0 -> BCMMaur1.0, mouse GRCm39 -> GRCm38
##############################################################################

# Input files:
#   /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_peaks.bed
#   /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_summits.bed
#
# liftOver binary and flags used below:
#   /Xnfs/lbmcdb/Semon_team/jcaron/utils/liftOver
#   -preserveInput
#   -bedPlus=3
#
# General liftOver usage: liftOver oldFile map.chain newFile unMapped

# Ham
/Xnfs/lbmcdb/Semon_team/jcaron/utils/liftOver \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_peaks_no_colnames.bed \
/Xnfs/lbmcdb/Semon_team/jcaron/files/chains/MesAur2.0_to_1.0.chain \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_peaks_no_colnames_lift_to_maur1.bed \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_peaks_no_colnames_lift_unmapped.bed \
-preserveInput \
-bedPlus=3

/Xnfs/lbmcdb/Semon_team/jcaron/utils/liftOver \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_summits_no_colnames.bed \
/Xnfs/lbmcdb/Semon_team/jcaron/files/chains/MesAur2.0_to_1.0.chain \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_summits_no_colnames_lift_to_maur1.bed \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_genebank_names_summits_no_colnames_lift_unmapped.bed \
-preserveInput \
-bedPlus=3

# Mus
/Xnfs/lbmcdb/Semon_team/jcaron/utils/liftOver \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_peaks_no_colnames_no_chr.bed \
/Xnfs/lbmcdb/Semon_team/jcaron/files/chains/GRCm39_to_GRCm38.chain \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_peaks_no_colnames_no_chr_lift_to_mm38.bed \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_peaks_no_colnames_no_chr_lift_unmapped.bed \
-preserveInput \
-bedPlus=3

/Xnfs/lbmcdb/Semon_team/jcaron/utils/liftOver \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_summits_no_colnames_no_chr.bed \
/Xnfs/lbmcdb/Semon_team/jcaron/files/chains/GRCm39_to_GRCm38.chain \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_summits_no_colnames_no_chr_lift_to_mm38.bed \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/mus_by_molar_with_best_summit_median_with_id_summits_no_colnames_no_chr_lift_unmapped.bed \
-preserveInput \
-bedPlus=3

##############################################################################
# 6. halLiftover step
##############################################################################

sbatch halliftover_job.sh

##############################################################################
# 7. Move the peak/summit name column into column 4 (standard BED name field)
##############################################################################

awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, $3, $6, $4, $5}' /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_peak_name_col4.bed

awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, $3, $6, $4, $5}' /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4.bed

awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, $3, $7, $4, $5, $6}' /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_summits_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_summits_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4.bed

##############################################################################
# 8. Run HALPER's orthologFind.py to call the final ortholog set
##############################################################################

python /Xnfs/lbmcdb/Semon_team/jcaron/utils/halLiftover-postprocessing/orthologFind.py \
-max_len 1000 -min_len 50 \
-protect_dist 5 \
-qFile /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_peak_name_col4.bed \
-tFile /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4.bed \
-sFile /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_summits_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4.bed \
-oFile /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered.bed \
-mult_keepone \
-preserve signal

##############################################################################
# 9. Reformat HALPER output and intersect with the mouse peak set
##############################################################################

awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, $3, $5, $4, $6, $7, $8, $9}' ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered.bed \
> ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered_bed_shape.bed

sed -E 's/^chr([0-9]+)\t/\1\t/' ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered_bed_shape.bed > ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered_bed_shape_no_chr.bed

# Overlap with >= 20% reciprocal coverage
bedtools intersect -a ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered_bed_shape_no_chr.bed \
-b mus_by_molar_with_best_summit_median_with_id_peaks_no_colnames_no_chr_lift_to_mm38.bed \
-f 0.2 \
-r \
> intersect.bed

# Overlap with any (>= 1bp) overlap
bedtools intersect -a ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered_bed_shape_no_chr.bed \
-b mus_by_molar_with_best_summit_median_with_id_peaks_no_colnames_no_chr_lift_to_mm38.bed \
> intersect_1bp.bed

wc -l intersect.bed

##############################################################################
# 10. Annotate peaks with their classification and tally counts per category
##############################################################################

awk 'BEGIN{FS=OFS="\t"} 
NR==FNR {val[$4]=$6; next} 
{print $0, ($4 in val ? val[$4] : "NA")}' \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_peak_name_col4.bed \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered_bed_shape_no_chr.bed > /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered_bed_shape_no_chr_classif.bed

cut -f10 /Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_lift_to_mm38_peak_name_col4_halpered_bed_shape_no_chr_classif.bed | sort | uniq -c

cut -f6 ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_peak_name_col4.bed | sort | uniq -c
  
awk 'BEGIN{FS=OFS="\t"} 
NR==FNR {val[$4]=$6; next} 
{print $0, ($4 in val ? val[$4] : "NA")}' \
/Xnfs/abc/nf_scratch/mrouyer/halper/files/ham_by_molar_with_best_summit_median_with_id_peaks_no_colnames_lift_nflo_chainfile_to_maur1_peak_name_col4.bed \
intersect.bed > intersect_classif.bed

cat intersect_classif.bed | cut -f10 | sort | uniq -c

cut -f5 mus_by_molar_with_best_summit_median_with_id_peaks_no_colnames_no_chr_lift_to_mm38.bed | sort | uniq -c

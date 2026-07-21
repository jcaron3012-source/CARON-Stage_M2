#!/bin/bash

#Cleans up a raw BED file for further analyses

bedfile=$1

echo "Removing column 5 from $bedfile"
cut -f 1-4,6- $bedfile > "${bedfile}.tmp" && mv "${bedfile}.tmp" "$bedfile"

echo "Replacing 'Distal Intergenic' with 'Distal_Intergenic'"
sed -i 's/Distal Intergenic/Distal_Intergenic/g' $bedfile

echo "Replacing ', intron X of X' with '_intron_X_of_X'"
sed -i -E 's/, (intron|exon) ([0-9]{1,2}) of ([0-9]{1,2})/_\1_\2_of_\3/g' $bedfile

echo "Replacing '3' UTR' with '3_UTR'"
sed -i "s/3' UTR/3_UTR/g" $bedfile

echo "Replacing '5' UTR' with '5_UTR'"
sed -i "s/5' UTR/5_UTR/g" $bedfile

echo "All replacements done."

#!/bin/bash

#converts a BED file to the right genome version, uses a chainfile.

bedfile=$1
chainfile=$2

#Main command
echo "liftOver of "$bedfile""
/Xnfs/lbmcdb/Semon_team/jcaron/utils/liftOver \
  --bedPlus \
  $bedfile \
  $chainfile \
  /Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/"liftOver_${bedfile}" \
  /Xnfs/lbmcdb/Semon_team/jcaron/files/unMapped/unMapped_links.bed
echo ""

echo "Smoothing up notations ..."
#Adding a chr prefix for further analyses
sed -i 's/^X/chrX/' "liftOver_${bedfile}"
sed -i 's/^\([0-9]\{1,2\}\)\s/chr\1 /' "liftOver_${bedfile}"
sed -i 's/^Y/chrY/' "liftOver_${bedfile}"

awk '/^chr/' "liftOver_${bedfile}" > "mus_${bedfile}"

echo "All done"

#!/bin/bash

#Extracts a specific genome from a HAL file, and converts it to a fasta

HALFILE=$1
SPECIES=$2
OUTFILE=$3

/home/cactus/bin/hal2fasta $HALFILE $SPECIES > $OUTFILE

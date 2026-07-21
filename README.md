### jcaron01 - CIGOGNE Team

------------
# Pipeline

![](Images/ATAC_ZOONOMIA_PIPELINE.png)

------------
# Repository structure

. \
├── Images \
│   ├── ATAC_ZOONOMIA_PIPELINE.png \
│   └── RodentTree.png \
├── mus_peaks.bed \
├── README.md \
├── repoTree.txt \
├── scripts \
│   ├── pipelines \
│   │   ├── conv.sh \
│   │   ├── hal2fasta.sh \
│   │   ├── hal2maf.sh \
│   │   ├── halExtract.sh \
│   │   ├── halLiftOver.sh \
│   │   ├── halStats.sh \
│   │   ├── msa_view_bulk.sh \
│   │   ├── msa_view.sh \
│   │   ├── phyloFit.sh \
│   │   └── RERtree.sh \
│   ├── python_scripts \
│   │   ├── check_seq_lengths.py \
│   │   ├── merge_msa_by_region_bulk.py \
│   │   └── merge_msa_by_region.py \
│   └── slurmjobs \
│       ├── slurm_conv.sh \
│       ├── slurm_hal2fasta.sh \
│       ├── slurm_hal2maf.sh \
│       ├── slurm_halExtract.sh \
│       ├── slurm_halLiftOver.sh \
│       ├── slurm_halStats.sh \
│       ├── slurm_msa_view_bulk.sh \
│       ├── slurm_msa_view.sh \
│       ├── slurm_phyloFit_command.sh \
│       ├── slurm_phyloFit.sh \
│       ├── slurm_phyloFit_summary.sh \
│       └── slurm_RERtree.sh \
└── trees \
    ├── ancien \
    │   └── tree.txt \
    ├── LowEpi100_tree.txt \
    ├── LowMes100_tree.txt \
    ├── PleioEpi100_tree.txt \
    ├── PleioMes100_tree.txt \
    ├── UpEpi100_tree.txt \
    └── UpMes100_tree.txt

8 directories, 37 files


------------

# Rodent tree

![](Images/helder-tree.png)

------------
# Main steps

 ## Step 1: Downloading alignment data from Zoonomia

_In node s92node01, in Xnfs/abc/nf_scratch/jcaron_

 ```bash
wget https://cgl.gi.ucsc.edu/data/cactus/447-mammalian-2022v1.hal -O zoonomia.hal 
 ```

 **File created:**

 - zoonomia.hal

 ## Step 2: Filtering species

 _In Xnfs/lbmcdb/Semon_team/jcaron_

 Use the Docker image `quay.io/comparative-genomics-toolkit/cactus:v3.1.4` (https://github.com/ComparativeGenomicsToolkit/cactus)

 **Code used:** 
 - `scripts/halExtract.sh`
 - `slurmjobs/slurm_halExtract.sh`

 **File created:**

 - rodents.hal

 ## Step 3: Conversion to MAF format

  _In scratch/Lake_

  Uses the Docker image `quay.io/comparative-genomics-toolkit/cactus:v3.1.4`

 **Code used:** 
 - `scripts/hal2maf.sh`
 - `slurmjobs/conv.sh`

  **Prerequisites:**
  - HAL file from Step 2
  - BED file (Martin's work) converted to the correct version (see `conv.sh`)
  - To convert the file, use UCSC's `liftOver` (https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver) and the ENSEMBL chain file (https://ftp.ensembl.org/ pub/current_assembly_chain/mus_musculus/GRCm39_to_GRCm38.chain.gz)

  **File created:**

  - rodents.maf

  ## Step 4: Converting the different alignments to FASTA format
  _In scratch/_

   Uses the **PHAST** package (https://github.com/UCSantaCruzComputationalGenomicsLab/phast/tree/master), specifically the `msa_view` function (https://manpages.debian.org/testing/phast/msa_view.1.en.html)

  **Code used:** 
 - `scripts/msa_view.sh`
 - `slurmjobs/slurm_msa_view.sh`
 - `python_scripts/merge_msa_by_region.py`

  **Prerequisites:**

 -  MAF file from Step 3

 **Files created:**

  - All .fa files in `fasta_msa`
 
## Step 5: Model Fitting
  _In scratch/_

  Uses the **PHAST** package (https://github.com/UCSantaCruzComputationalGenomicsLab/phast/tree/master), specifically the `phyloFit` function (https://manpages.debian.org/testing/phast/phyloFit.1.en.html)

  **Code used:** 
 - `scripts/phyloFit.sh`
 - `slurmjobs/slurm_phyloFit.sh`
  **Prerequisites:**

 -  Fasta files from Step 4, in `fasta_msa`

 **Files created:**

  - All .mod files in `modfiles`

   ## Step 6: Extracting the trees 
  _In scratch/_

  **Code used:** 
 - `scripts/RERtree.sh`

  **Prerequisites:**

 -  Mod files from Step 5, in `modfiles`

 **File created:**

  - tree.txt in `trees`

## Step 7: RERconverge analysis

_In RStudio_

**Code used:**
- `Rscripts/rerconverge.R`

---------
# Optional steps
 ## Extracting a specific genome   from a HAL file

 **Code used:**
 - `scripts/hal2fasta.sh`

 **Prerequisites**
 - `halfiles/rodents.hal`

 **File created**
 - `mm.fa`

 ## Accessing info from the HAL file

 **Code used:**
 - `scripts/halStats.sh`

 **Prerequisites**
 - `halfiles/rodents.hal`

 **Files created**
 - `halStats/`

  ## Plotting the phylogenetic tree and alignments

 **Code used:**
 - `Rscripts/plot_phylo.R`

 **Prerequisites**
 - `plot_phylo/rodents.nh`
 - `plot_phylo/HelderPhenotypes.csv`

 **Files created**
 - `Images/helder-tree.png`
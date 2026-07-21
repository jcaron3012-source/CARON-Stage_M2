#!/bin/bash
#SBATCH -J hal2maf
#SBATCH -o /Xnfs/lbmcdb/Semon_team/jcaron/log/hal2maf_%j_%a.out
#SBATCH -e /Xnfs/lbmcdb/Semon_team/jcaron/log/hal2maf_%j_%a.err
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH -p Lake
#SBATCH -t 1-00:00:0

#Uses a WGA in HAL format and a BED file, aligns all the sequences present in the BED file on the WGA.

set -ex

echo -n "Time started: "
date

species_name=$1
hal_file_input=$2
output_dir=$3
chunk_dir=$4
remove_dups=$5

#Separate the BED file into chunks for parallelisation
chunk_id=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")
bed_chunk="${chunk_dir}/chunk_${chunk_id}.bed"

if [ ! -f "$bed_chunk" ]; then
	echo "Chunk file $bed_chunk not found"
	exit 1
fi

maf_file="${output_dir}/chunk_${chunk_id}.maf"

if [ -s "${maf_file}" ]; then
	echo "${maf_file} already exists and is not empty. Skipping."
	echo -n "Time ended: "
	date
	exit 0
fi

WORKDIR="/scratch/Lake/jcaron/tmp/work_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
toil_home="/scratch/Lake/jcaron/tmp/toil_home_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
jobstore="/scratch/Lake/jcaron/jobstore_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORKDIR" "$toil_home"

CHRUN=/Xnfs/abc/charliecloud_bin/ch-run
IMG=/Xnfs/abc/charliecloud/img/cactus_v3.1.4.sqfs

options="--refGenome ${species_name} --onlyOrthologs --noAncestors --outType single \
--bedRanges ${bed_chunk} --binariesMode local --workDir ${WORKDIR} --batchCores 4 --maxDisk 400G"

#Main command: cactus-hal2maf on every chunk
$CHRUN \
	--write \
	--bind=/scratch/:/scratch/ \
	--bind=/Xnfs:/Xnfs \
	--bind=$toil_home:/root/.toil \
	$IMG \
	-- bash -c "export PATH=/home/cactus/cactus_env/bin:/home/cactus/bin:\$PATH && cactus-hal2maf $jobstore $hal_file_input $maf_file $options"

rm -rf "$jobstore" "$toil_home" "$WORKDIR"

#Remove the duplicates if wanted
if [[ "$remove_dups" == 1 ]]; then
	/Xnfs/lbmcdb/Semon_team/jcaron/utils/mafTools/bin/mafDuplicateFilter -m "$maf_file" > "${maf_file}_tmp"
	mv "${maf_file}_tmp" "$maf_file"
fi

gzip "$maf_file_output"

echo -n "Time ended: "
date

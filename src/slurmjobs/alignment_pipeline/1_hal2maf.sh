#!/bin/bash
#SBATCH -J submit-hal2maf
#SBATCH -o /Xnfs/lbmcdb/Semon_team/jcaron/log/submit_hal2maf_%j.out
#SBATCH -e /Xnfs/lbmcdb/Semon_team/jcaron/log/submit_hal2maf_%j.err
#SBATCH -n 1
#SBATCH -t 00:10:00
set -ex

#Submit the hal2maf job.

echo -n "Time started: "
date

species_name="Mus_musculus"
hal_file_input="/Xnfs/lbmcdb/Semon_team/jcaron/files/halfiles/rodents_noamericans.hal"
coordinates_file="/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/mus_links.bed"
output_dir="/Xnfs/lbmcdb/Semon_team/jcaron/files/maffiles"
chunk_dir="/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/chunks"
error_dir="/Xnfs/lbmcdb/Semon_team/jcaron/log"
array_job_script="/Xnfs/lbmcdb/Semon_team/jcaron/src/pipelines/alignment_pipeline/1-src_hal2maf_v2.sh"
remove_dups=1 

#Parallelization parameters
desired_jobs=200  
throttle=40      

if [ ! -f "$hal_file_input" ]; then
	echo "No HAL file"
	exit 1
fi
if [ ! -f "$coordinates_file" ]; then
	echo "No coordinates file"
	exit 1
fi

mkdir -p "$output_dir" "$chunk_dir" "$error_dir" || exit 1

#Checking # chunks
existing_chunks=$(ls "$chunk_dir"/chunk_*.bed 2>/dev/null | wc -l)

#Determining # jobs
if [ "$existing_chunks" -eq 0 ]; then
	total_regions=$(wc -l < "$coordinates_file")
	max_array_size=$(scontrol show config 2>/dev/null | awk -F= '/MaxArraySize/{gsub(/ /,"",$2); print $2}')
	max_array_size=${max_array_size:-1000}
	max_array_size=$(( max_array_size - 1 ))

	num_chunks=$desired_jobs
	if [ "$num_chunks" -gt "$max_array_size" ]; then
		num_chunks=$max_array_size
	fi
	if [ "$num_chunks" -gt "$total_regions" ]; then
		num_chunks=$total_regions
	fi
	lines_per_chunk=$(( (total_regions + num_chunks - 1) / num_chunks ))

	echo "Total regions: $total_regions | MaxArraySize: $((max_array_size+1)) | chunks: $num_chunks | lines/chunk: $lines_per_chunk"

	awk -v n="$lines_per_chunk" -v dir="$chunk_dir" '
	{
		idx = int((NR-1)/n) + 1
		printf "%s\n", $0 > (dir "/chunk_" sprintf("%04d", idx) ".bed")
	}
	' "$coordinates_file"
else
	echo "Chunks already in $chunk_dir ($existing_chunks files)"
fi

num_chunks_actual=$(ls "$chunk_dir"/chunk_*.bed 2>/dev/null | wc -l)
if [ "$num_chunks_actual" -eq 0 ]; then
	echo "Chunk splitting failed, no chunk files produced."
	exit 1
fi

#Listing jobs to do
missing=()
for f in "$chunk_dir"/chunk_*.bed; do
	base=$(basename "$f" .bed)         
	id=${base#chunk_}                   
	maf="${output_dir}/${base}.maf"
	if [ ! -s "$maf" ]; then
		missing+=("$((10#$id))")
	fi
done

if [ "${#missing[@]}" -eq 0 ]; then
	echo "All is already converted."
	echo -n "Time ended: "
	date
	exit 0
fi

IFS=$'\n' sorted=($(sort -n <<<"${missing[*]}")); unset IFS

array_spec=""
range_start=${sorted[0]}
prev=${sorted[0]}
for ((i=1; i<${#sorted[@]}; i++)); do
	cur=${sorted[$i]}
	if [ "$cur" -eq $((prev+1)) ]; then
		prev=$cur
		continue
	fi
	if [ "$range_start" -eq "$prev" ]; then
		array_spec="${array_spec}${range_start},"
	else
		array_spec="${array_spec}${range_start}-${prev},"
	fi
	range_start=$cur
	prev=$cur
done
if [ "$range_start" -eq "$prev" ]; then
	array_spec="${array_spec}${range_start}"
else
	array_spec="${array_spec}${range_start}-${prev}"
fi

echo "Chunks missing (${#missing[@]}/${num_chunks_actual}): $array_spec"

if [ "$throttle" -gt 0 ]; then
	array_spec="${array_spec}%${throttle}"
fi

#Launch script
sbatch --array="$array_spec" -p Lake "$array_job_script" \
	"$species_name" "$hal_file_input" "$output_dir" "$chunk_dir" "$remove_dups"

echo -n "Time ended: "
date

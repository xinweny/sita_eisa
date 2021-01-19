#!/bin/bash
filename=$1

while IFS="" read -r GSE || [ -n "$GSE" ]
do
    echo "submitting job for $GSE"
    command="sbatch --export=GSE=$GSE,threads=$threads \
    --job-name=MAPCOUNT-$GSE \
    --cpus-per-task=$threads \
    02QuasRAlignmentAndCounting.sbatch"
    echo $command
    $command
done < "$filename"

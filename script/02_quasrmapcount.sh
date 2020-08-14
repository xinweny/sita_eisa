#!/bin/bash
filename=$1

while IFS="" read -r GSE || [ -n "$GSE" ]
do
    echo "submitting job for $GSE"
    command="sbatch --export=GSE=$GSE --job-name=$GSE-MAPCOUNT -o ./log/$GSE_aligncount.out 02QuasRAlignmentAndCounting.sbatch"
    echo $command
    $command
done < "$filename"

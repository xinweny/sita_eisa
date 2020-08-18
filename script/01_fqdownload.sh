#!/bin/bash
filename=$1

while IFS="" read -r GSE || [ -n "$GSE" ]
do
    echo "submitting job for $GSE"
    command="sbatch --export=GSE=$GSE \
    --job-name=FQ-$GSE \
    01FastqDownload.sbatch"
    echo $command
    $command
done < "$filename"

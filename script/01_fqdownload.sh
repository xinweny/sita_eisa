#!/bin/bash
filename=$1

while IFS="" read -r GSE || [ -n "$GSE" ]
do
    echo "submitting job for $GSE"
    command="sbatch --export=GSE=$GSE --job-name=$GSE-FQDOWNLOAD -o ./log/$GSE_fqdownload.out 01FastqDownload.sbatch"
    echo $command
    $command
done < "$filename"

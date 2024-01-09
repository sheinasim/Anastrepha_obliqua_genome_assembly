#!/bin/bash

echo "Starting assembly at: $(date)"

version=`hifiasm --version`

echo "Using HiFiASM version: ${version}"

/project/ag100pest/sheina.sim/software/hifiasm/hifiasm -o ${prefix} -t 48 --primary *fcsfilt.fastq.gz

echo "Finished assembling at: $(date). Printing stats."

for x in `ls *_ctg.gfa | sed 's/\.gfa//g'`; do echo ${x}; any2fasta ${x}.gfa > ${x}.fasta; done

for x in `ls *fasta | sed 's/\.fasta//'`
do
stats.sh -Xmx4g ${x}.fasta >${x}.stats
done

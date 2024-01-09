#!/bin/bash

mkdir tmp

ls *fcsfilt.fastq.gz > FILES

kmc -k21 -t48 -m64 -ci1 -cs10000000 @FILES reads tmp/

kmc_tools transform reads histogram kmcreads10000000.histo -cx10000000

genomescope2 -i kmcreads10000000.histo -n ${species} -o ${outdir} -p ${ploidy} 

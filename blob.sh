#!/bin/bash

workingdir=`pwd | tr "/" "\t" | awk '{print $(NF)}'`

if stat --printf='' *csi 2>/dev/null
then
	echo "Incomplete blobdir exists, deleting blobdir so it can be made again."
	rm -rf ${assembly}
	rm *csi
fi

if [[ -d "std" ]]
then
   echo "std directory exists."
else
   mkdir std
fi

if [[ ! -s "${assembly}.reads.bam" ]]
then
	echo "Mapping reads to assembly."
	minimap2 -ax map-hifi -t 96 ${assembly}.fasta all_files.fcsfilt.fastq.gz | samtools sort -@96 -O BAM -o ${assembly}.reads.bam -
else 
	echo "Reads are mapped to the assembly."
fi

if [[ ! -s "${assembly}_blast.out" ]]
then
	echo "Performing nucleotide blast of assembly."
	blastn -db $nucleotide_db/nt \
        -query ${assembly}.fasta \
        -outfmt "6 qseqid staxids bitscore std" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads 48 \
        -out ${assembly}_blast.out
else
	echo "Nucleotide blast finished."
fi

if [[ ! -s "${assembly}_records_shorter_than_10000000.fasta" ]]
then
	echo "Creating .fasta with records shorter than 10mb."
	python /project/ag100pest/sheina.sim/software/seqer/short_seqer.py 10000000 ${assembly}.fasta
fi

if [[ ! -s "${assembly}_diamond.out" ]]
then
	echo "Performing uniprot blast of assembly using diamond blast."
	diamond blastx \
        --query ${assembly}_records_shorter_than_10000000.fasta \
        --db $blob_dbs/uniprot/reference_proteomes.dmnd \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads 48 \
        > ${assembly}_diamond.out
else
	echo "Uniprot blast complete."
fi

buscodir=run_`ls ${assembly}_auto_euk | grep "specific" | sed 's/short_summary.specific.//g' | sed "s/${assembly}.*//" | sed 's/odb10\./odb10/'`

busco_table=${assembly}_auto_euk/${buscodir}/full_table.tsv

if [[ ! -d "${assembly}" ]]
then
	echo "Creating blob dir."
    blobtools create \
        --fasta ${assembly}.fasta \
        ${assembly}

    blobtools add \
        --cov ${assembly}.reads.bam \
        ${assembly}

    blobtools add \
        --hits ${assembly}_blast.out \
		--hits ${assembly}_diamond.out \
		--taxrule bestsumorder \
		--taxdump ${blob_dbs}/taxdump \
		${assembly}
    
	blobtools add --busco ${busco_table} ${assembly}
else
	echo "Blob directory complete."
fi

if [[ ! -s "${assembly}_blobblurbout.tsv" ]]
then
	echo "Creating blobblurb summary."
    python $blurb_path/blobblurb.py ${assembly} ${busco_table} ${assembly}.fasta
else
	echo "Blobblurb summary exists."
fi

## Create figures

if [[ ! -d "${assembly}_figures" ]]
then
	mkdir ${assembly}_figures
fi

sed -i '/level/d' ${assembly}/meta.json

singularity exec $softwarepath/Singularity_files/blobtk_latest.sif blobtk plot -v snail -d ${assembly} -o ${assembly}_figures/${assembly}_scaff_snail.png
singularity exec $softwarepath/Singularity_files/blobtk_latest.sif blobtk plot -v blob -d ${assembly} -o ${assembly}_figures/${assembly}_scaff_blob.png
singularity exec $softwarepath/Singularity_files/blobtk_latest.sif blobtk plot -v cumulative -d ${assembly} -o ${assembly}_figures/${assembly}_scaff_cumulative.png

## Identify contaminant records

cat ${assembly}_blobblurbout.tsv | grep -v "Arthropoda\|record" | awk -v OFS='\t' '{print $1}' >${assembly}_blob_contaminant_contigs.txt

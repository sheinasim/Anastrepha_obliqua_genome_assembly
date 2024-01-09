#!/bin/bash

busco \
	-i ${assembly}.fasta \
	-o ${assembly}_auto_euk \
	--auto-lineage-euk \
	-m geno \
	-c 48 

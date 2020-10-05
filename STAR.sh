#!/bin/bash

STAR="/data/mg1/caix/src/biosoft/STAR-2.6.1d/bin/Linux_x86_64/STAR"

mkdir -p star_index

${STAR}    --runMode genomeGenerate \
           --runThreadN 10  \
           --genomeDir star_index \
           --genomeFastaFiles /data/mg1/caix/Pangenome_analysis/CCB_rna/genome.fa \
        #--sjdbGTFfile /home/share/genome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf \
        #--sjdbOverhang 149

${STAR}    --runThreadN 20 \
           --runMode alignReads \
           --genomeDir star_index \
           --readFilesIn HN53_1_1.fq.ft.gz,HN53_2_1.fq.ft.gz,HN53_3_1.fq.ft.gz,HN53_4_1.fq.ft.gz HN53_1_2.fq.ft.gz,HN53_2_2.fq.ft.gz,HN53_3_2.fq.ft.gz,HN53_4_2.fq.ft.gz   \
           --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outWigType wiggle read2















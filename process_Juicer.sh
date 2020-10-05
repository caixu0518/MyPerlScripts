#!/bin/bash 

export PATH=/home/caix/miniconda3/envs/py2.7/bin/:${PATH}

sample=JZS_V2.180805
draftgenome=Bol-JZS.PacBio.final_180623.fasta.anchoredScf.fa
chrsize=Bol-JZS.PacBio.final_180623.fasta.anchoredScf.fa.sizes

#bwa index $draftgenome
#python  /data/mg1/caix/src/Hi-C/Juicer/juicer-1.6.2/misc/generate_site_positions.py HindIII $sample  $draftgenome

#then, put all raw fastq file to fastq folder (mkdir fastq) and then make sure you have changed the fastq file name(i.e. *R1.fastq.gz *R2.fastq.gz)
/data/mg1/caix/src/Hi-C/Juicer/juicer-1.6.2/CPU/juicer.sh  -g $sample -s HindIII -z $draftgenome -y $sample"_HindIII.txt" -p  $chrsize -t 40 -D /data/mg1/caix/src/Hi-C/Juicer/juicer-1.6.2/CPU  &> juicer.log

#!/bin/bash

sample="Bol"
draftgenome="Bol-jzs.PacBio.scf.fa.corrected_scf.fa"
chrsize="Bol-jzs.PacBio.scf.fa.corrected_scf.fa.sizes"

##------##process Juicer

#bwa index ${draftgenome}
#python /data/mg1/caix/src/Hi-C/Juicer/juicer-1.6.2/misc/generate_site_positions.py  HindIII ${sample}  ${draftgenome}

#then, put all raw fastq file to fastq folder (mkdir fastq) and then make sure you have changed the fastq file name(i.e. *R1.fastq.gz *R2.fastq.gz)
bash  /data/mg1/caix/src/Hi-C/Juicer/juicer-1.6.2/CPU/juicer.sh  -g ${sample} -s HindIII -z ${draftgenome} -y ${sample}"_HindIII.txt" -p  ${chrsize} -t 80 -D  /data/mg1/caix/src/Hi-C/Juicer/juicer-1.6.2/CPU  &> juicer.log

##------##process 3d-dna

#path=`pwd`

#output=${sample}".3d-dna"

#mkdir ${output}
#cd ${output}

#ln -sf ${path}"/aligned/merged_nodups.txt" .
#ln -sf ${path}"/"${draftgenome} .

#bash /data/mg1/caix/src/Hi-C/3d-dna/run-asm-pipeline.sh -m haploid -r 0  ${draftgenome}  "merged_nodups.txt"  &> 3d_dna.log

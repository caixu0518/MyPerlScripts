#!/bin/bash

sample="B.nigra"
read1="B.nigra-250bp_1.fq.gz"
read2="B.nigra-250bp_2.fq.gz"

########################
readLen=150
kmer=21
threads=20
hashsize="25G"
genomescope="/data/mg1/caix/src/biosoft/genomescope/genomescope.R"

temp1=${sample}".Temp_1.fq"
temp2=${sample}".Temp_2.fq"
gzip -dc ${read1} > ${temp1} &&  gzip -dc ${read2} > ${temp2} 

jellyfish count -C -m ${kmer}  -s ${hashsize}  -c 8  -t  ${threads} ${temp1} ${temp2}  -o ${sample}".reads.jf"

jellyfish histo -t ${threads}  ${sample}".reads.jf" -o  ${sample}".reads.histo"

Rscript ${genomescope} ${sample}".reads.histo" ${kmer} ${readLen} ${sample}".output" &> ${sample}".kmer.Stat"

rm ${temp1} ${temp2} ${sample}".reads.jf"

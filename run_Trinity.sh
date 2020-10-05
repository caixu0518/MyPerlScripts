#!/bin/bash -ve

export PATH=/home/caix/miniconda3/envs/py2.7/bin/:${PATH}
#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

#/data/mg1/caix/src/biosoft/trinityrnaseq-Trinity-v2.6.6/Trinity --seqType fq --max_memory 1000G \
#              --samples_file sample.txt \
#              --CPU 44


/data/mg1/caix/src/biosoft/trinityrnaseq-Trinity-v2.6.6/Trinity  \
                --genome_guided_bam rnaseq_alignments.csorted.bam \
                --max_memory 100G  \
                --genome_guided_max_intron 10000 --CPU 20 \

##### Done Running Trinity #####

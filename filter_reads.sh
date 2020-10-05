#!/bin/bash

##---sample list
SAMPLES=sample.txt

##-begin to filter reads
samples=$(cat $SAMPLES);

for sample in ${samples}
do

fastp -i ${sample}"_1.fq.gz"  -I ${sample}"_2.fq.gz"  -o ${sample}"_1.fq.ft.gz" -O ${sample}"_2.fq.ft.gz"  -z 4 -q 20 -u 30 -n 5 -w 10 -h ${samples}.html

done

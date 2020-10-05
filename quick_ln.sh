#!/bin/bash

SAMPLES=/data/mg1/caix/works/Bjuncea_Pop_analysis/splitGenomes/samples.list
samples=$(cat $SAMPLES)

for sample in ${samples}
do 
ln -s /data/mg1/caix/works/Bjuncea_Pop_analysis/Bju_Pop_clean_reads/${sample}* .
done

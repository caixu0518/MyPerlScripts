#!/bin/bash 

genome1=""
genome2=""

###---###---
folder=`pwd`
scanPAV_version="/data/mg1/caix/src/biosoft/scanPAV/src"

cpus=30
aligner=bwa
score=550  
output=$aligner\_$score


# Assembly locations:
assembly1=$folder/$genome1
assembly2=$folder/$genome2
outname=$(basename $assembly1 .fasta)


rm -rf $output
mkdir -p $output 
cd $output

ln -sf $assembly1 .
ln -sf $assembly2 .

echo; echo Looking for Presence PAVs | tee log.txt ; echo
$scanPAV_version/scanPAV -nodes $cpus -score $score -align $aligner $(basename $assembly1) $(basename $assembly2)  pavs_present_in_$outname | tee -a  log.txt
echo; echo Looking for Absence PAVs | tee -a log.txt ; echo
$scanPAV_version/scanPAV -nodes $cpus -score $score -align $aligner $(basename $assembly2) $(basename $assembly1)  pavs_absent_in_$outname  | tee -a  log.txt

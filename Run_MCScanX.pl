#!/usr/bin/perl -w
use strict;
##- Xu Cai

my $in0 = $ARGV[0];     ##- query gene coords
my $in1 = $ARGV[1];     ##- reference gene coords
my $in2 = $ARGV[2];     ##- query pep
my $in3 = $ARGV[3];     ##- reference pep 
my $prefix = $ARGV[4];  ##- prefix e.g. Z1_to_Tair10

my $MCScanX = '/data/mg1/caix/src/biosoft/MCScanX';

##- change coords to MCScanX format
my $querygff = "query.gff";
   `awk '{print \$2"\t"\$1"\t"\$3"\t"\$4}' $in0 > $querygff`;

my $refgff = "ref.gff";
   `awk '{print \$2"\t"\$1"\t"\$3"\t"\$4}' $in1 > $refgff`;

##- merge gff files
system("cat $querygff $refgff > $prefix.gff");

my $allprotein = 'all.proteins.fasta';
   `cat $in2 $in3 > $allprotein`;

##-----Run blastp--------
system("makeblastdb -in $allprotein  -parse_seqids -hash_index  -out $allprotein  -dbtype prot");
system("blastp -query $allprotein  -db $allprotein  -out  $prefix.blast -outfmt 6  -num_threads 18 -num_alignments 5   -evalue 1e-10");
print "Step 1: blasp finished\n";

##-----Run MCScanX-------
system("mkdir MCScanX_data");

system("mv  $prefix.gff  $prefix.blast  MCScanX_data");

print "Step 2: begin to run MCScanX\n";
system("$MCScanX/MCScanX MCScanX_data/$prefix");

system("$MCScanX/duplicate_gene_classifier MCScanX_data/$prefix");

##- clean 
system("rm $querygff $refgff  $allprotein $allprotein.*");
print "all finished\n";

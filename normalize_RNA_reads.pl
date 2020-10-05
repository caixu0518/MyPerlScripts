#!/usr/bin/perl -w
use strict;

my ($left, $right) = ($ARGV[0], $ARGV[1]);

my $wd="/data/mg1/caix/src/biosoft/trinityrnaseq-Trinity-v2.6.6/util";
my $normalize_outdir = "insilico_read_normalization";

## Normalize reads
my $max_memory = "20G";
my $max_read_coverage = 50;
my $min_kmer_cov = 1;
my $normalize_max_pct_stdev = 10000;


my $cmd = "$wd/insilico_read_normalization.pl --seqType fq --JM $max_memory "
        . " --max_cov $max_read_coverage --min_cov $min_kmer_cov --CPU $CPU --output $normalize_outdir   --max_pct_stdev $normalize_max_pct_stdev "
        . "--left $left --right $right";








































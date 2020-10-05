#!/usr/bin/perl -w
##- Xu Cai
use strict;

my $subreads = $ARGV[0]; ##- merged.raw.subreads.bam
my $primers = $ARGV[1]; ##- primers file

my $cpus = 20;
my $cmdString = 'NA';

##- begin CCS ...
    $cmdString = "ccs $subreads  ccs.bam  --minSnr 2   --noPolish --minPasses 1  --numThreads $cpus";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

##- begin lima ...
    $cmdString = "lima  --isoseq --dump-clips --no-pbi -j $cpus  ccs.bam  $primers  demux.bam";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

##- begin cluster ...
    $cmdString = "isoseq3 cluster demux.primer_5p--primer_3p.bam  unpolished.bam  --verbose --require-polya  -j $cpus";    
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

##- begin polish ...
    $cmdString = "isoseq3 polish -j $cpus  unpolished.bam  $subreads   polished.bam";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

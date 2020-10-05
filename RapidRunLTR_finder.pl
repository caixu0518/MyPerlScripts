#!/usr/bin/perl -w
use strict;
use threads;

my $in0 = $ARGV[0]; #genome file

my %id2seq = ();
open IN0, $in0;
while(<IN0>){
chomp;
 if(/^>(\S+)/ && ! exists $id2seq{$1}){
  $id2seq{$1} = "";
 }
 else{
  $id2seq{$1} .= $_;
 }
}
close IN0;

my @filename = ();
open OUT0, ">scaffold.fa";
for my $key1(sort keys %id2seq){
 print $key1, "\n";
 if($key1 =~ /scaffold/){   #it depends__
   print OUT0 ">$key1","\n", $id2seq{$key1}, "\n";
 }
 else{
  push(@filename, "$key1.fa");
  open OUT1, ">$key1.fa";
  print OUT1 ">$key1", "\n", $id2seq{$key1}, "\n";
  close OUT1;
 }
}
close OUT0;
push(@filename, "scaffold.fa");

my @thr = ();
for(my $i = 0; $i <= $#filename; $i+= 1){  
  $thr[$i] = threads->create(\&LTR_finder, $filename[$i]);
} 

for(my $i = 0; $i <= $#filename; $i+= 1){
  $thr[$i]->join;
}


sub LTR_finder{
 my($input) = @_;
 system("ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9  -s /share/mg6t/caix/src/LTR_FINDER.x86_64-1.0.6/tRNAdb/Athal-tRNAs.fa  $input > $input.ltr_finder");
} 

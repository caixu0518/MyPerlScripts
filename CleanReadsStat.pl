#!/usr/bin/perl -w
use strict;

open IN, $ARGV[0]; #compressed reads file (Pair-end reads name were seperated by Tab)
open OUT, ">ReadsStat.list";
open OUT5, ">Unmatched.txt";
while(<IN>){
chomp;
 my @a = split("\t", $_);
 if(-f $a[0] && -f $a[1]){
    my @b = &readsStat($a[0], $a[1]);
    print OUT join("\t", @b),"\n";
 }
 else{
    print OUT5 join("\t", @a),"\n";
 }
}
close IN;
close OUT;


sub readsStat {
  my ($L, $R) = @_;
  $L =~ /^(\S+?)_1\.fq\.ft\.gz$/;  #it depends__
  my ($flag, $StatL, $StatR, $readsData, $readL) = (0, 0, 0, 0, 0);
  open (IN0, "gzip -dc $L |");
  while(<IN0>){
  chomp;
  $flag += 1;
    $StatL += length($_), if($flag == 2);
    $flag =0, if($flag == 4);
  }
  close IN0;
 
 open (IN1, "gzip -dc $R |");
 while(<IN1>){
 chomp;
 $flag += 1;
   $StatR += length($_), if($flag == 2);
   $flag = 0, if($flag == 4);
   $readL = length($_), if(eof);
 }
 close IN1;
 $readsData = $StatL + $StatR;
 my @result = ("$1", $readL, $readsData);
 return(@result);
}

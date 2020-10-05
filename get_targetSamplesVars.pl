#!/usr/bin/perl -w
#Xu Cai
use strict;

my $in0 = $ARGV[0]; ##- sample
my $in1 = $ARGV[1]; ##- position
my $in2 = $ARGV[2]; ##- genotype
my $in3 = $ARGV[3]; ##- target sample list


    &main();
sub main {

  my %targetSamples = ();
     &read_targetSamples($in3, \%targetSamples);

  my $out0 = $in0.".root";
  my $out1 = $in1.".root";
  my $out2 = $in2.".root";
     &output($in0, $in1, $in2, \%targetSamples, $out0, $out1, $out2);     

}

sub output {

   my ($sample, $position, $genotype, $targetSamples, $sampleft, $positionft, $genotypeft) = @_;

   my @file = ();
   my $count = 0;
   open IN1, $sample;
   open OUT1, ">$sampleft";
   while(<IN1>){
     chomp;
     if(exists $targetSamples ->{$_}){
        print OUT1 $targetSamples ->{$_}, "\n";
        push(@file, $count);
        
     }
     $count += 1;    
   }
   close IN1;
   close OUT1;

   ##- process genotype
   open IN2, $genotype;
   open OUT2, ">$genotypeft";
   while(<IN2>){
     chomp;
     $_ =~ s/([AGCTN]{2})/$1\t/g;
     my @a = split("\t", $_);
     my @b = ();
     for my $key(@file){
         push(@b, $a[$key]);
     }
     print  OUT2 join("", @b),"\n";
   }
   close IN2;
   close OUT2;

   ##- process position file
   `cp $position $positionft`;

}

sub read_targetSamples {
 
    my ($targetFile, $targetSamples) = @_;

    open IN0, $targetFile;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
      if(not exists $targetSamples ->{$temp[0]}){
         $targetSamples ->{$temp[0]} = $temp[1];
      }
      else{
         print "Warinings: same id $temp[0]\n";
      }
    }
    close IN0;

}

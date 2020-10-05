#!/usr/bin/perl -w
##- Xu Cai
use strict;

my $in0 = $ARGV[0]; ##- chromosome sizes
my $in1 = $ARGV[1]; ##- coordinate file(bed format)
my $index = $ARGV[2]; ##- index 
my $window = 500000;
my $step = 50000;

    &main();
sub main {

   my %chr2Len = ();
      &read_chr($in0, \%chr2Len);
  
   my $out = $index."_windows.list";
      &generate_value(\%chr2Len, $in1, $window, $step, $out);


}


sub generate_value {

    my ($chr2Len, $bedFile, $windowSize, $stepSize, $out) = @_;
    
    my $indexSNP;
    open IN0, $bedFile;
    while(<IN0>){
      chomp;
      my ($chr, $pos, $value) = (split(/\t/, $_))[0,1,3];
          $indexSNP ->{$chr} ->{$pos} = $value;
    }
    close IN0;
   
    open OUT, ">$out";        
    for my $chrid(sort keys %{$chr2Len}){
        my $num_step = 0;
        my ($num, $all, $flag) = (0,0,0);
        
        for($num_step = 0; ($num_step*$step + $window) <= $chr2Len ->{$chrid}; $num_step++){
            ($num, $all) = (0, 0);
            my ($start, $end) = ($num_step*$step, $num_step*$step + $window -1);
            for my $pos(sort {$a<=>$b} keys %{$indexSNP ->{$chrid}}){
                if($pos >= $start && $pos <= $end){
                   $all += 1;
                   $num += 1, if($indexSNP ->{$chrid} ->{$pos} == 1);
                }
                last, if($pos > $end);
            }
            $flag = 0, if($all == 0);
            $flag = $num/$all, if($all > 0);
            print OUT $chrid, "\t", $num_step*$step,  "\t", ($num_step*$step + $window -1), "\t" , $flag, "\n";         
        }         
        if($num_step*$step < $chr2Len ->{$chrid} && ($num_step*$step + $window) > $chr2Len ->{$chrid}){
           ($num, $all) = (0, 0);
           my ($start, $end) = ($num_step*$step, $chr2Len ->{$chrid});
           for my $pos(sort {$a<=>$b} keys %{$indexSNP ->{$chrid}}){
               if($pos >= $start && $pos <= $end){ 
                  $all += 1; 
                  $num += 1, if($indexSNP ->{$chrid} ->{$pos} == 1);
               }
           } 
           $flag = 0, if($all == 0);
           $flag = $num/$all, if($all > 0);
           print OUT  $chrid, "\t", $num_step*$step,  "\t", ($num_step*$step + $window), "\t" , $flag, "\n";
        }
    }
    close OUT;

}

sub read_chr {

    my ($chrSizes, $chr2Len) = @_;

    open IN0, $chrSizes;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
      $chr2Len ->{$temp[0]} = $temp[1];
    }
    close IN0;

}


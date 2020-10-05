#!/usr/bin/perl -w
use strict;
use threads;

my $in0 = $ARGV[0]; #sample file

my %sam2index = ();
my $sampleN = 0;
   $sampleN = &sampleindex($in0, \%sam2index);

my $out = "$in0.Stat";
   &output(\%sam2index, $out);

###---###---###---###
sub output {
  my ($samindex, $out) =@_;
  my (@thr);

  for(my $i=1; $i<=$sampleN; $i+=1) {
    $thr[$i] = threads->create(\&ReadsStat, $samindex ->{$i});
  }
  for(my $i=1; $i<=$sampleN; $i+=1) {
    $thr[$i]->join;
  }
  
  system("rm -rf $out") if(-f $out);
  for(my $m=1; $m<=$sampleN; $m+=1){
    my $statF = $samindex ->{$m}.".Stat";
    system("cat $statF >> $out");
    system("rm -rf $statF");
  }
                                                                                                                                        
}

sub sampleindex {
  my ($samF, $samindex) = @_;
  my $n = 0;

  open IN0, $samF;
  while(<IN0>){
    $n += 1;
    if(/(\S+)/){
      $samindex ->{$n} = $1;
    }
  }
  close IN0;
  return($n);
}


sub ReadsStat {
  my ($sam) = @_;
  my ($L) = $sam."_1.fq.ft.gz"; #it depends_
  my ($R) = $sam."_2.fq.ft.gz"; #it depends_
  my @Lresult = &singleStat($L);
  my @Rresult = &singleStat($R);
  my @Stat = ($sam, $Lresult[0] + $Rresult[0], $Lresult[1]);   #sample, reads stat, read length
 
  open OUT0, ">$sam.Stat";
  print OUT0 join("\t", @Stat), "\n";
  close OUT0;
}

sub singleStat {
  my ($ReadF) = @_;
  my ($stat, $flag, $readlen) = (0, 0, 0);  

  open (IN0, "gzip -dc $ReadF |");
  while(<IN0>){
    chomp;
    $flag += 1;
    $stat += length($_), if($flag == 2);
    $flag =0, if($flag == 4);
    $readlen = length($_), if(eof);
  }
  close IN0;
  my @result = ($stat, $readlen);
  return(@result);
}

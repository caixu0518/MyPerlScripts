#!/usr/bin/perl -w
#Author: Feng Cheng
use warnings;
use strict;

my $in = $ARGV[0];
my $col = $ARGV[1];
my $start = $ARGV[2];
my $inter = $ARGV[3];
my $out = $in.".hist.col".$col.".step".$inter;

&process($in,$col,$start,$inter,$out);

sub process {
  my ($in,$col,$start,$inter,$out) = @_;

  my @data;
  &readData($in,\@data,$col);

  my ($min,$max);
  my @datahist = &histData(\@data,$inter,$start,\$min,\$max);

  &output($out,\@datahist,$max);

}

sub output {
  my ($out,$datahist,$maxy) = @_;

  my $maxx = -100;
  my @int2freq;
  open(my $FW,">$out");
  print $FW "Start\tStop\tNumber\tFrequence\n";
  for(my $i=0;$i<@$datahist;$i+=1) {
    print $FW $datahist->[$i],"\n";
    my @data = split(/\t/,$datahist->[$i]);
    if($data[2]>$maxx) {
      $maxx = $data[2];
    }
  }

  print $FW "\n\n";
  for(my $i=0;$i<@$datahist;$i+=1) {
    my @data = split(/\t/,$datahist->[$i]);
    $data[0] = sprintf "%0.2f", $data[0];
    $data[1] = sprintf "%0.2f", $data[1];

    my $int = "$data[0]-$data[1]";
    my $lenInt = length($int);
    my $lenSpa = " "x(10-$lenInt);
    print $FW $int,$lenSpa;
    for(my $j=1;$j<=int($data[3]*1000+0.5);$j+=1) {
      print $FW "*";
    }
    if($data[2]==$maxx) {
      print $FW " --> $data[2]";
    }
    print $FW "\n";
  }
  close($FW);

}

sub readData {
  my ($in,$data,$col) = @_;

  open(my $FR,$in);
  my $index = 0;
  while($_=<$FR>) {
    my @temp = split(/\t/);
    if($temp[$col]=~/^[-+\d\.]+$/) {
      #if($temp[$col] < 10) {
        $data->[$index] = $temp[$col];
        $index += 1;
      #}
    }
  }
  #my @data = map{(split(/\t/))[2]} <FR>;
  close($FR);

}


sub histData {
  my $nargv = @_;
  my ($data,$interval,@histed,$minStart,$min,$max);
  if($nargv==5) {
    ($data,$interval,$minStart,$min,$max) = @_;
  }
  elsif($nargv==4) {
    ($data,$interval,$min,$max) = @_;
  }
  $$min = &min($data);
  #$minimum = 0;   # it depends

  my $size = @$data;
  for(my $i=0;$i<$size;$i+=1) {
    if($data->[$i] eq "NA") {
      $data->[$i] = -1;
    }
  }

  $$max = &max($data);
  print "$minStart\t$interval\t$$min\t$$max\n";
  my $index = 0;
  for(my $i=$$min+$interval;($i-$interval)<=$$max;$i+=$interval) {
    my $start = $i-$interval;
    my $num = &find_num_between($data,$start,$i);
    my $freq = $num/$size;
    $histed[$index] = $start."\t".$i."\t".$num."\t".$freq;
    $index += 1;
  }
  return @histed;
}

sub max {
  my ($data) = @_;
  my $max = $data->[0];
  for(my $i=1;$i<@$data;$i+=1) {
    if($max<$data->[$i]) {
      $max = $data->[$i];
    }
  }
  return $max;
}

sub min {
  my ($data) = @_;
  my $min = $data->[0];
  for(my $i=1;$i<@$data;$i+=1) {
    if($min>$data->[$i]) {
      $min = $data->[$i];
    }
  }
  return $min;
}

sub find_num_between {
  my ($data,$start,$stop) = @_;
  my $num = 0;
  for(my $i=0;$i<@$data;$i+=1) {
    if($data->[$i]>=$start && $data->[$i]<$stop) {
      $num += 1;
    }
  }
  return $num;
}

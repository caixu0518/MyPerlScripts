#!/usr/bin/perl -w

#Program: Feng Cheng
#Bug report: chengfeng@caas.cn

use strict;
use warnings;
use threads;
use Getopt::Long;



=head1 program discription

        Call Variations From Population Samples With High Accuracy

        Prerequiried tools for running: BWA, Samtools


=head1 Command-line Option

        --inDir              [required] The directory where the population resequencing reads are stored.
        --outDir             [required] The directory for store the filtered reads and other temp file.
        --qual               [required] File to provide the illumina reads quality
                                          version 1.8 (all >=1.8), or 1.3 (all <1.8).
        --ref                [required] Fasta format of the reference genome.
        
	--outStat            [optional] File to store the statistics during filtering raw reads.
        --Gsize              [optional] The size of the genome in units of Gb. 
                                          Automatic calculated from the reference genome as default.
        --folds              [optional] How many folds for each sample use to do the calling.
                                          Use all data as default [0].
        --build              [optional] Build index [1] for the reference, or don't build it [0] 
                                          incase it has been built previously. [1].
        --pop                [optional] The file of polymorphic loci for variant calling.
                                          If not exist, will perform pooled mapping as default [1] 
        --sam                [optional] File to store the samples used in the variant calling.
        --pos                [optional] File to store the genomic position of variants called.
        --geno               [optional] File to store the genotypes of variants called.
        --help


=head1 Usage

        perl popSeq2Geno.pl -inDir example_reads/ -outDir example_reads_temp/ -qual example_quality -ref example_genome.fa -folds 1 -build 1 -pop 1

=cut

my $softpath = "/home/caix/bin";

##--Initialize input data--##

my ($dir,$outDir,$qual,$stat,$gsize,$fold,$refG,$bref,$pop,$sam,$pos,$geno,$help);

GetOptions(
        "inDir=s"=>\$dir,
        "outDir=s"=>\$outDir,
        "qual=s"=>\$qual,
        "outStat=s"=>\$stat,
        "Gsize=s"=>\$gsize,
        "folds=s"=>\$fold,
        "ref=s"=>\$refG,
        "build=s"=>\$bref,
        "pop=s"=>\$pop,
        "sam=s"=>\$sam,
        "pos=s"=>\$pos,
        "geno=s"=>\$geno,
        "help"=>\$help,
);


die `pod2text $0` if($help);
die `pod2text $0` if(!defined $dir);
system("rm -rf Reads_temp") if(-e "Reads_temp");
if(! defined $outDir){
 $outDir = "Reads_temp/";
}
die `pod2text $0` if(!defined $outDir);
die `pod2text $0` if(!defined $qual);
die `pod2text $0` if(!defined $refG);


if(!defined $stat) {
  $stat = "readsFilter.statistics";
}

if(!defined $gsize) {
  $gsize = &readGsize($refG);
}

if(!defined $fold) {
  $fold = 0;
}

if(!defined $bref) {
  $bref = 1;
}

if(!defined $pop) {
  $pop = 1;
}

if(!defined $sam) {
  $sam = "samples";
}
if(!defined $pos) {
  $pos = "positions.loci";
}
if(!defined $geno) {
  $geno = "genotypes.loci";
}

#if($pop !~ /^1$/) {
#  $fold = 0;
#}


#--Launch the running process-----------------------

&process_main($dir,$outDir,$qual,$stat,$gsize,$fold,$refG,$bref,$pop,$sam,$pos,$geno);


sub process_main {

  my ($dir,$outDir,$qual,$stat,$gsize,$fold,$ref,$bref,$pop,$sam,$pos,$geno) = @_;
  
  my (%file,%rfile2rec);
  &readFile(\%file,\%rfile2rec,$dir);

  if(!(-e $outDir)) {
    system("mkdir $outDir");
  }

  my $fOutDir = $outDir."reads_ft/";
  &qualityFilter($qual,\%file,\%rfile2rec,$dir,$fOutDir,$stat);

  print "\nAll filtering process run successfully...\n\n";


  my %sam2fold;
  &readStat($stat,\%sam2fold,$gsize);

  my $pfOutDir = $outDir."reads_sel/";
  &pickOutReads($fOutDir,$pfOutDir,$fold,\%sam2fold);
  

  if($bref == 1) {
    system("$softpath/bwa index $ref");       #--format reference--
    system("$softpath/samtools faidx $ref");
  }
  
  my $mapDir = $outDir."mapped/";
  my $plist;
  if($pop =~ /^1$/) {
    $plist = "populationPolymorphicLoci.posAlleleT";

    &bMapPop($pfOutDir,$mapDir,$ref,$plist,$fold,\%sam2fold);
  }
  else {
    $plist = $pop;
  }

  my $mapDir2 = $mapDir;
  if($fold == 0) {
    &bMapIndi($pfOutDir,$mapDir2,$ref,$plist,$pos,$geno,$sam,$fold,\%sam2fold,$pop);
  }
  else {
    $mapDir2 = $outDir."mapped2/";
    &bMapIndi($fOutDir,$mapDir2,$ref,$plist,$pos,$geno,$sam,$fold,\%sam2fold,$pop);
  }

}


sub readStat {

  my ($stat,$sam2fold,$gsize) = @_;

  open(my $FR,$stat);
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);

    $temp[0] =~ /\/([^\/]+)_1\.fq[s]*\.gz/;
    my $sam = $1;
    my $rfold = $temp[7]/1000000000/$gsize;
    
    $sam2fold->{$sam} = $rfold;
    #print "$sam\t$rfold\n";
  }
  close($FR);

}


sub pickOutReads {

  my ($fOutDir,$pfOutDir,$fold,$sam2fold) = @_;

  if(-e $pfOutDir) {
    system("rm -rf $pfOutDir");
  }
  system("mkdir $pfOutDir");

  if($fold!=0) {
    foreach my $sam (sort {$a cmp $b} keys %$sam2fold) {
      my $rr1 = $sam."_1.fq.ft.gz";
      my $rr2 = $sam."_2.fq.ft.gz";
      my $r1 = $sam."_1.fq.gz";
      my $r2 = $sam."_2.fq.gz";
      my $in1 = $fOutDir.$rr1;         #caix
      my $in2 = $fOutDir.$rr2;         #caix
      my $out1 = $pfOutDir.$r1;      #caix
      my $out2 = $pfOutDir.$r2;      #caix
      if($sam2fold->{$sam}>$fold) {
        my $nfold = $fold/$sam2fold->{$sam};
        
        open(my $FR1,"gzip -dc $in1 |");
        open(my $FR2,"gzip -dc $in2 |");
        open(my $FW1,"| gzip >$out1");
        open(my $FW2,"| gzip >$out2");

        my (@L1,@L2);
        while($_=<$FR1>) {
          $L1[0] = $_;
          $L1[1] = <$FR1>;
          $L1[2] = <$FR1>;
          $L1[3] = <$FR1>;
          
          
          $L2[0] = <$FR2>;
          $L2[1] = <$FR2>;
          $L2[2] = <$FR2>;
          $L2[3] = <$FR2>;
          if(rand(1)<$nfold) {
            print $FW1 $L1[0];
            print $FW1 $L1[1];
            print $FW1 $L1[2];
            print $FW1 $L1[3];

            print $FW2 $L2[0];
            print $FW2 $L2[1];
            print $FW2 $L2[2];
            print $FW2 $L2[3];
          }
        }

        close($FR1);
        close($FR2);
        close($FW1);
        close($FW2);
      }
      else {
        system("ln -s ../../$in1 $out1");
        system("ln -s ../../$in2 $out2");  #caix
      }
    }
  }
  else {
    foreach my $sam (sort {$a cmp $b} keys %$sam2fold) {
      #my $r1 = $sam."_1.fq.ft.gz";
      #my $r2 = $sam."_2.fq.ft.gz";
      #my $ls1 = $fOutDir.$r1;
      #my $ls2 = $fOutDir.$r2;
      #$ls1 =~s/^[^\/]+\///;
      #$ls2 =~s/^[^\/]+\///;
      #system("ln -s ../$ls1 $pfOutDir");
      #system("ln -s ../$ls2 $pfOutDir");
      my $rr1 = $sam."_1.fq.ft.gz";
      my $rr2 = $sam."_2.fq.ft.gz";
      my $r1 = $sam."_1.fq.gz";
      my $r2 = $sam."_2.fq.gz";
      my $in1 = $fOutDir.$rr1;  #caix
      my $in2 = $fOutDir.$rr2;  #caix
      my $out1 = $pfOutDir.$r1;
      my $out2 = $pfOutDir.$r2;
      system("ln -s ../../$in1 $out1");
      system("ln -s ../../$in2 $out2"); #caix
    }
  }
}


sub readGsize {

  my ($in) = @_;

  my %sca2seq;
  &readFasta($in,\%sca2seq);

  my $gsize = 0;
  foreach my $sca (keys %sca2seq) {
    $gsize += length($sca2seq{$sca});
  }
  $gsize /= 1000000000;

  return $gsize;

}


sub readFasta {

  my ($in,$id2seq) = @_;
  open(my $SFR,$in);

  my $id;
  while($_=<$SFR>) {
    if(/^>([^\s^\n]+)\s*\n*/) {
      $id = $1;
      $id2seq->{$id} = "";
    }
    else {
      chomp;
      $id2seq->{$id} .= $_;
    }
  }
  close($SFR);

}


sub bMapPop {

  my ($dir,$outDir,$ref,$plist,$fold,$sam2fold) = @_;

  my (%file,%file2rec);
  ##print "$dir\n";
  &readFile(\%file,\%file2rec,$dir);

  &bwaMapPop(\%file,\%file2rec,$dir,$outDir,$ref,$plist,$fold,$sam2fold);

  print "\nPopulation mapping process run successfully...\n\n";

  &mergeBam($outDir,$plist,$ref,$fold,$sam2fold);

}


sub bwaMapPop {

  my ($file,$file2rec,$dir,$outDir,$ref,$plist,$fold,$sam2fold) = @_;

  if(-e $outDir) {
    system("rm -rf $outDir");
  }
  system("mkdir $outDir");

  my $sample;
  my $posF = $plist;  #--it depends----

  foreach $sample (sort {$a cmp $b} keys %$file) {

    print "\nStart processing $sample...\n";
   # my $in1 = $sample."_1.fq.ft.gz";
   # my $in2 = $sample."_2.fq.ft.gz";
    
    my $in1 = $sample."_1.fq.gz";   #caix
    my $in2 = $sample."_2.fq.gz";    
  
    #print "$in1\t$in2\n";
    if(exists($file2rec->{$in1}) && exists($file2rec->{$in2})) {

      my $outSam = $outDir.$sample.".cf.sam";

      &bmap($dir,$outDir,$in1,$in2,$outSam,$sample,$ref);

      system("rm $outSam");
    }
  }

}


sub mergeBam {

  my ($outDir,$plist,$ref,$fold,$sam2fold) = @_;

  my $outDirm = $outDir."*.bam";

  my $samID = "merge_mapped";
  system("$softpath/samtools merge -@ 6  $samID.bam $outDirm");

  my $t1 = $samID.".bam.sorted";
  system("$softpath/samtools sort  $samID.bam $t1");
  #system("rm $samID.bam");
  print "bam sorting is OK!\n";

  my $t2 = $t1.".bam";
  system("$softpath/samtools index $t2");

  my @sam = keys %$sam2fold;
  my $nsam = @sam;

  my $bdep;
  if($fold==0) {
    foreach my $sam (keys %$sam2fold) {
      $bdep += $sam2fold->{$sam};
    }
  }
  else {
    foreach my $sam (keys %$sam2fold) {
      if($sam2fold->{$sam}<$fold) {
        $bdep += $sam2fold->{$sam};
      }
      else {
        $bdep += $fold;
      }
    }
  }

  my ($minDep,$maxDep);
  #$minDep = $bdep * 0.5;
  #if($minDep < 10) {
  #  $minDep = 10;
  #}

  $maxDep = $bdep * 2;
  $maxDep = $bdep * 4;
  #print "$maxDep: $ref: $t2: $samID\n";
  system("$softpath/samtools mpileup -q 30 -Q 30 -C50 -d $maxDep -L $maxDep -ugf $ref $t2 | $softpath/bcftools view -p 0.9999 -vcg - > $samID.cons.vcf");

  my $vfile = $samID.".cons.vcf";
  #print "$vfile\n";

  my $minAdep = 3;         #--it depends---
  &screenGenotypeBwa($vfile,$maxDep,$minAdep);

  $minDep = 10;             #--it depends---
  system("awk '\(\$5!~/,/ && \$5~/\\w/ && \$6+\$7>=$minDep\)' $samID.cons.vcf.genotype > $plist");

  #system("rm $t2 $t2.bai");
  #system("rm $vfile $samID.cons.vcf.genotype");

}


sub screenGenotypeBwa {

  my ($in,$dep,$minAdep) = @_;
  my $out = $in.".genotype";

  open(my $FR,$in);
  open(my $FW,">$out");

  ##min depth for alleles in total samples is set here!

  while($_=<$FR>) {

    chomp;
    my @temp = split(/\t/);
    #if($temp[4]=~/\w/ && !($temp[4]=~/,/)) {
    #if(/^[^#]/ && $temp[4]=~/\w/) {
    if(/^[^#]/) {
      if($temp[7] =~ /;DP4=(\d+),(\d+),(\d+),(\d+);/) {
        my ($d1,$d2) = ($1+$2,$3+$4);
        my ($ref,$a1,$a2) = ($temp[3],$temp[3],$temp[4]);

        if($d1>=$minAdep && $d1<=$dep && $d2>=$minAdep && $d2<=$dep) {
          #if($a2 eq $ref) {
          #  $a2 = $a1;
          #  $a1 = $ref;
          #}
          #elsif($a1 ne $ref && $a2 ne $ref) {
          #  $a1 = "N";
          #  $a2 = "N";
          #}
        }
        elsif($d1>=$minAdep && $d1<=$dep && $d2==0) {
          $a2 = $a1;
          $d2 = $d1;
        }
        elsif($d2>=$minAdep && $d2<=$dep && $d1==0) {
          $a1 = $a2;
          $d1 = $d2;
        }
        else {
          $a1 = "N";
          $a2 = "N";
        }
        #if(!($a2=~/^\w$/)) {
        #  $a1 = "N";
        #  $a2 = "N";
        #}

        if($a1 ne "N" && $a2 ne "N" && $a1 ne $a2) {
          print $FW "$temp[0]\t$temp[1]\t$ref\t$a1\t$a2\t$d1\t$d2\n";
        }

      }
    }
  }

  close($FR);
  close($FW);

}



sub bMapIndi {

  my ($dir,$outDir,$ref,$plist,$posiF,$genoF,$samF,$fold,$sam2fold,$pop) = @_;

  #if($pop !~ /^1$/) {
  #  if(-e $outDir) {
  #    system("rm -rf $outDir");
  #  }
  #  system("mkdir $outDir");
  #}

  my (%file,%file2rec);
  &readFile(\%file,\%file2rec,$dir);

  &bwaMap(\%file,\%file2rec,$dir,$outDir,$ref,$plist,$posiF,$genoF,$samF,$fold,$sam2fold,$pop);

  print "\nAll individual calling process run successfully...\n\n\n";

}


sub bwaMap {

  my ($file,$file2rec,$dir,$outDir,$ref,$plist,$posiF,$genoF,$samF,$fold,$sam2fold,$pop) = @_;

  if(!(-e $outDir)) {
    system("mkdir $outDir");
  }

  ##my ($sample,%i2geno);
  my $sample;
  my $posF = $plist;  #--it depends----

  ##-snp--
  ##my (%i2pos,%i2alle);
  ##&readPos($posF,\%i2pos,\%i2alle);

  foreach $sample (sort {$a cmp $b} keys %$file) {

    print "\nStart processing $sample...\n";
    my ($in1,$in2);
    #$in1 = $sample."_1.fq.ft.gz";
    #$in2 = $sample."_2.fq.ft.gz";
    $in1 = $sample."_1.fq.gz";
    $in2 = $sample."_2.fq.gz";

    if(exists($file2rec->{$in1}) && exists($file2rec->{$in2})) {
      my $outSam = $outDir.$sample.".cf.sam";

      if($fold!=0 || $pop!~/^1$/) {
        &bmap($dir,$outDir,$in1,$in2,$outSam,$sample,$ref);
      }

      ##&callVar($outSam,$ref,\%i2geno,$posF,\%i2pos,\%i2alle);
      &callVar($outSam,$ref,$posF);

    }
  }

  ##system("perl extractGP.pl $outDir $posF $posiF $genoF $samF");

  print "\nWe are here for outputing genotypes...\n";
  &extractGP($outDir,$posF,$posiF,$genoF,$samF,$fold,$sam2fold,);

  ##print "$outDir\n";

}


sub extractGP {
  
  my ($outDir,$posF,$posiF,$genoF,$samF,$fold,$sam2fold) = @_;

  #-snp--
  my (%i2pos,%i2alle);
  &readPos($posF,\%i2pos,\%i2alle);

  my %i2geno;
  &extractGeno($outDir,\%i2pos,\%i2alle,\%i2geno,$samF,$sam2fold);

  my ($posiF1,$genoF1) = ($posiF.".snp",$genoF.".snp");
  &output(\%i2pos,\%i2alle,\%i2geno,$posiF1,$genoF1);

  #-indel--
  my (%pindex2oldPos,%oldPos2rec,%ppos2indel);
  &popIndel($posF,\%pindex2oldPos,\%oldPos2rec,\%ppos2indel);

  my %pos2indelG;
  &extractIndel($outDir,\%pindex2oldPos,\%oldPos2rec,\%ppos2indel,\%pos2indelG,$samF,$sam2fold);

  my ($posiF2,$genoF2) = ($posiF.".indel",$genoF.".indel");
  &outputIndel(\%pindex2oldPos,\%oldPos2rec,\%ppos2indel,\%pos2indelG,$posiF2,$genoF2);

}


sub outputIndel {
  my ($index2oldpos,$oldPos2rec,$ppos2indel,$pos2geno,$posiF,$genoF) = @_;

  open(my $FW1,">$posiF");
  open(my $FW2,">$genoF");
  
  my %pt2rec;
  foreach my $i (sort {$a <=> $b} keys %$index2oldpos) {
    my $opos = $index2oldpos->{$i};
    my $pt = $oldPos2rec->{$opos};
    if(!exists($pt2rec{$pt})) {
      $pt2rec{$pt} = "Y";
      if(exists($pos2geno->{$pt})) {
        my $geno = $pos2geno->{$pt};
        my $tg = $geno;
        $tg =~ s/N//g;
        if($tg=~/\w/) {
          print $FW1 "$ppos2indel->{$pt}\n";
          print $FW2 "$geno\n";
        }
      }
    }
  }
  close($FW1);
  close($FW2);

}


sub extractIndel {

  my ($outDir,$pindex2oldPos,$oldPos2rec,$ppos2indel,$pos2geno,$samF,$sam2fold) = @_;

  open(my $FW,">$samF.indel");
  opendir(DIR,$outDir) || die "Error in opening dir $outDir\n";
  my $filename;

  my $mdep;
  while(($filename = readdir(DIR))){
    if($filename=~/^(.+)\.cf\.vcf\.list/) {
      my $sam = $1;
      print $FW "$sam\n";

      my $rfold = $sam2fold->{$sam};
      if($rfold>3) {
        $mdep = 0;                        #**to be revised***************
      }
      else {
        $mdep = 0;
      }

      my %pos2rec;
      &readIndel($outDir.$filename,\%pos2rec);

      my %indel2rec;
      my %crtPt2geno = ();
      my $index = 1;
      foreach my $index (keys %$pindex2oldPos) {
        my $opos = $pindex2oldPos->{$index};
        my $pt = $oldPos2rec->{$opos};
        my $indel = $ppos2indel->{$pt};
        if(!exists($crtPt2geno{$pt}) || $crtPt2geno{$pt} ne "GG") {
          if(exists($pos2rec{$opos})) {
            my $data = $pos2rec{$opos};
            if($data =~ /INDEL/) {
              my @temp = split(/\t/,$data);
              my ($d1,$d2);
              if($temp[7]=~/;DP4=(\d+),(\d+),(\d+),(\d+);/) {
                ($d1,$d2) = ($1+$2,$3+$4);
              }

              my ($spt,$sindel);
              my ($len1,$len2) = (length($temp[3]),length($temp[4]));
              &typeIndel($temp[0],$temp[1],$temp[3],$temp[4],$len1,$len2,\$spt,\$sindel);
              if($indel eq $sindel) {
                if($d1>0 && $d2>0) {
                  $crtPt2geno{$spt} = "CG";
                }
                elsif($d2>0) {
                  $crtPt2geno{$spt} = "GG";
                }
                elsif($d1>0) {
                  $crtPt2geno{$spt} = "CC";
                }
                else {
                  if(!exists($crtPt2geno{$spt})) {
                    $crtPt2geno{$spt} = "NN";
                  }
                }
              }
              else {
                if(!exists($crtPt2geno{$pt})) {
                  $crtPt2geno{$pt} = "NN";
                }
              }
            }
            else {
              $crtPt2geno{$pt} = "CC";
            }
          }
          else {
            if(!exists($crtPt2geno{$pt})) {
              $crtPt2geno{$pt} = "NN";
            }
          }
        }
      }

      my %pt2rec = ();
      foreach my $index (keys %$pindex2oldPos) {
        my $opos = $pindex2oldPos->{$index};
        my $pt = $oldPos2rec->{$opos};
        if(!exists($pt2rec{$pt})) {
          $pt2rec{$pt} = "Y";
          my $cgeno = $crtPt2geno{$pt};
          if(exists($pos2geno->{$pt})) {
            $pos2geno->{$pt} .= $cgeno;
          }
          else {
            $pos2geno->{$pt} = $cgeno;
          }
        }
      }
    }
  }
  close($FW);
}


sub readIndel {
  my ($in,$pos2rec) = @_;

  open(my $FR,$in);
  while($_=<$FR>) {
    if(!(/^#/)) {
      chomp;
      my @temp = split(/\t/);

      if(!($temp[4]=~/,/)) {
        my $pos = $temp[0]."\t".$temp[1];
        $pos2rec->{$pos} = $_;
      }
    }
  }
  close($FR);
}


sub popIndel {
  my ($in,$pindex2oldPos,$oldPos2rec,$pos2indel) = @_;

  my $count = 1;
  my %pt2rec;
  open(my $FR,$in);
  open(my $FW,">$in.InDel");
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);

    my ($pold,$pt,$indel);
    my ($len1,$len2) = (length($temp[3]),length($temp[4]));
    if($len1>1 || $len2>1) {
      &typeIndel($temp[0],$temp[1],$temp[3],$temp[4],$len1,$len2,\$pt,\$indel);
      $pold = "$temp[0]\t$temp[1]";
      $pindex2oldPos->{$count} = $pold;
      $oldPos2rec->{$pold} = $pt;
      #$index2pos->{$count} = $pt;
      $pos2indel->{$pt} = $indel;
      $count += 1;
      if(!exists($pt2rec{$pt})) {
        print $FW "$indel\t$temp[5]\t$temp[6]\n";
        $pt2rec{$pt} = "Y";
      }
    }
  }
  #foreach my $i (sort {$a<=>$b} keys %$pindex2oldPos) {
  #  print "$i\t$pindex2oldPos->{$i}\t$oldPos2rec->{$pindex2oldPos->{$i}}\n";
  #}
  close($FR);
  close($FW);
}


sub typeIndel {
  my ($chr,$pos,$a1,$a2,$len1,$len2,$pt,$indel) = @_;

  my ($type,$size,$length,$allel);
  $length = abs($len1-$len2);
  $size = $length;
  if($len2<$len1) {
    $type = "D";
    $pos += $len2;
    $a1 =~ s/^(\w){$len2}//;
    $allel = $a1;
  }
  else {
    $type = "I";
    $pos += $len1;
    $a2 =~ s/^(\w){$len1}//;
    $allel = $a2;
  }
  $$pt = "$chr\t$pos";
  $$indel = "$$pt\t$type$size\t$allel";
}


sub output {

  my ($i2pos,$i2alle,$i2geno,$posiF,$genoF) = @_;

  my ($out1,$out2) = ($posiF,$genoF);   #--it depends--
  open(my $FW1,">$out1");
  open(my $FW2,">$out2");
  my $count = 0;

  foreach my $i (sort {$a <=> $b} keys %$i2pos) {
    my $alle = $i2alle->{$i};
    my @ta = split(/:/,$alle);
    my $tg = $i2geno->{$i};
    $tg =~ s/N//g;

    if($tg=~/\w/) {
      $tg =~ s/$ta[0]//g;
      $tg =~ s/$ta[1]//g;
      if(!($tg=~/\w/)) {
        my $pos = $i2pos->{$i};
        $pos =~ s/-/\t/;
        print $FW1 "$pos\t$ta[0]\t$ta[1]\n";
        print $FW2 "$i2geno->{$i}\n";
      }
    }
    else {
      $count += 1;
    }

  }

  #print "$count\n";
  close($FW1);
  close($FW2);

}


sub extractGeno {

  my ($outDir,$i2pos,$i2alle,$i2geno,$samF,$sam2fold) = @_;
  ##my $gfile;  #.list;

  open(my $FW,">$samF.snp");
  opendir(DIR,$outDir) || die "Error in opening dir $outDir\n";
  my $filename;

  my $mdep;
  while(($filename = readdir(DIR))){
    if($filename=~/^(.+)\.cf\.vcf\.list/) {
      my $sam = $1;
      print $FW "$sam\n";

      my $rfold = $sam2fold->{$sam};
      if($rfold>3) {
        $mdep = 0;                        #**to be revised***************
      }
      else {
        $mdep = 0;
      }

      my %pos2rec;
      &readData($outDir.$filename,\%pos2rec);
      foreach my $i (sort {$a <=> $b} keys %$i2pos) {
        my $pos = $i2pos->{$i};
        my @sa;
        ($sa[0],$sa[1]) = ("N","N");

        if(!exists($pos2rec{$pos})) {
          if(!exists($i2geno->{$i})) {
            $i2geno->{$i} = "NN";
          }
          else {
            $i2geno->{$i} .= "NN";
          }
        }
        else {
          my $data = $pos2rec{$pos};
          my @temp = split(/\t/,$data);
          if($temp[7]=~/;DP4=(\d+),(\d+),(\d+),(\d+);/) {
            my ($d1,$d2) = ($1+$2,$3+$4);
            if($d1>$mdep && $d2>$mdep) {
              if($temp[3]=~/^\w$/ && $temp[4]=~/^\w$/) {
                $sa[0] = $temp[3];
                $sa[1] = $temp[4];
              }
            }
            elsif($d1>$mdep) {
              if($temp[3]=~/^\w$/ && $temp[4]=~/^.{1}$/) {
                $sa[0] = $temp[3];
                $sa[1] = $temp[3];
              }
            }
            elsif($d2>$mdep) {
              if($temp[3]=~/^\w$/ && $temp[4]=~/^\w$/) {
                $sa[0] = $temp[4];
                $sa[1] = $temp[4];
              }
            }
          }
          if(!exists($i2geno->{$i})) {
            $i2geno->{$i} = "$sa[0]$sa[1]";
          }
          else {
            $i2geno->{$i} .= "$sa[0]$sa[1]";
          }
        }

      }
    }
  }

  close($FW);
  closedir(DIR);

}


sub readData {

  my ($gfile,$pos2rec) = @_;

  open(my $FR,$gfile);
  my $index = 1;

  while($_=<$FR>) {
    if(!(/^#/) && !(/INDEL/)) {
      chomp;
      my @temp = split(/\t/);

      if($temp[3]=~/^\w$/ && $temp[4]=~/^[\.\w]$/) {
        my $pos = $temp[0]."-".$temp[1];
        $pos2rec->{$pos} = $_;
      }

    }
  }
  close($FR);

}


sub readPos {

  my ($posF,$i2pos,$i2alle) = @_;

  open(my $FR,$posF);
  open(my $FW,">$posF.SNP");
  my $index = 1;
  my %pt2rec;
  while($_=<$FR>) {
    if(!(/^#/)) {
      chomp;
      my @temp = split(/\t/);
      if($temp[3]=~/^\w$/ && $temp[4]=~/^\w$/) {

        $i2pos->{$index} = $temp[0]."-".$temp[1];
        $i2alle->{$index} = $temp[3].":".$temp[4];

        $index += 1;
        my $pt = "$temp[0]\t$temp[1]";
        if(!exists($pt2rec{$pt})) {
          print $FW "$temp[0]\t$temp[1]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\n";
          $pt2rec{$pt} = "Y";
        }
      }
    }
  }
  close($FR);
  close($FW);
}


sub callVar {

  my ($outSam,$ref,$posF) = @_;

  $outSam =~ /^(.+)\.sam$/;
  my $sample = $1;
  system("rm $outSam");

  my $t1 = $sample.".bam";
  my $t2 = $t1.".sorted";
  system("$softpath/samtools sort -@ 20 $t1 $t2");
  system("rm $t1");

  my $t3 = $t1.".sorted.bam";
  system("$softpath/samtools index $t3");

  my $t4 = $sample.".vcf";
  system("$softpath/samtools mpileup -q 20 -Q 30 -C50 -l $posF -ugf $ref $t3 | $softpath/bcftools view -p 0.99 -cg - > $t4.list");

  ##system("rm $t3 $t3.bai");

  print "$sample calling is clear now...\n";

}


sub readMdep {

  my ($in,$sample2dep) = @_;

  open(my $FR,$in);
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);
    $sample2dep->{$temp[0]} = $temp[1];
  }
  close($FR);

}


sub bmap {

  my ($dir,$outDir,$in1,$in2,$out3,$sample,$ref) = @_;

  my $out1 = $outDir.$in1.".sai";
  my $out2 = $outDir.$in2.".sai";
  $in1 = $dir.$in1;
  $in2 = $dir.$in2;

  open(my $TT,$in1);
  $_=<$TT>;$_=<$TT>;chomp;
  my $rlength = length($_);
  close($TT);

  if($rlength>=70) {
    system("$softpath/bwa mem -t 5  $ref $in1 $in2 > $out3");
  }
  else {
    system("$softpath/bwa aln -t 6 $ref $in1 > $out1");
    system("$softpath/bwa aln -t 6 $ref $in2 > $out2");
    system("$softpath/bwa sampe -a 1000 $ref $out1 $out2 $in1 $in2 > $out3");
    system("rm $out1 $out2");
  }
  $out3 =~ /^(.+)\.sam$/;
  my $sam = $1;

  my $t1 = $sam.".bam";
  system("$softpath/samtools view -bS -@ 6  $out3 > $t1");

  print "$sam mapping is clear now...\n";

}


#sub readFile {
#
#  my ($file,$file2rec,$dir) = @_;
#
#  opendir(DIR,$dir) || die "Error in opening dir $dir\n";
#  my $filename;
#
#  while(($filename = readdir(DIR))){
#    if($filename=~/^(.+)_\d\.fq\.ft\.gz/) {     #it depends--- _hb_1.fq.ft _hb_2.fq.ft
#      $file->{$1} = "Y";
#      $file2rec->{$filename} = "Y";
#    }
#  }
#
#  closedir(DIR);
#
#}


sub qualityFilter {

  my ($qual,$file,$file2rec,$dir,$outDir,$stat) = @_;

  if(!(-e $outDir)) {
    system("mkdir $outDir");
  }

  if(-e $stat) {
    system("rm $stat");
  }
  
  my $sample;

  foreach $sample (sort {$a cmp $b} keys %$file) {
    #print "$sample\t";
    my $in1 = $sample."_1.fq.gz";   #--it depends--
    my $in2 = $sample."_2.fq.gz";

    if(exists($file2rec->{$in1}) && exists($file2rec->{$in2})) {

      #print "submit to filtering...\n";
      my $out1 = $outDir.$sample."_1.fq.ft";
      my $out2 = $outDir.$sample."_2.fq.ft";

      $in1 = $dir.$in1;
      $in2 = $dir.$in2;

      &filterPE($in1,$in2,$qual,$out1,$out2,$stat);

    }
    else {
      print("PE data not available...\n");
    }
  }

}


sub filterPE {

  my ($in1,$in2,$in3,$out1,$out2,$out3) = @_;

  my ($trimQ,$delAQ,$minLen,$Npercent) = (20,25,40,0.01);  #--for merge map   #--it depends--

  &processFilter($in1,$in2,$in3,$out1,$out2,$out3,$trimQ,$delAQ,$minLen,$Npercent);

}


sub processFilter {

  my ($in1,$in2,$in3,$out1,$out2,$out3,$trimQ,$delAQ,$minLen,$Npercent) = @_;

  my %asc2ps;
  my $headLink;

  my %sam2pool;
  &readSam($in3,\%sam2pool);

  $in1 =~ /([^\/]+)_1.fq/;
  #print "$in1\t$1\n";

  #-check code for quality
  #Sanger #!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ
  %asc2ps = ('!'=>0,'"'=>1,'#'=>2,'$'=>3,'%'=>4,
             '&'=>5,'\''=>6,'('=>7,')'=>8,'*'=>9,
             '+'=>10,','=>11,'-'=>12,'.'=>13,'/'=>14,
             '0'=>15,'1'=>16,'2'=>17,'3'=>18,'4'=>19,
             '5'=>20,'6'=>21,'7'=>22,'8'=>23,'9'=>24,
             ':'=>25,';'=>25,'<'=>27,"="=>28,'>'=>29,
             '?'=>30,'@'=>31,'A'=>32,'B'=>33,'C'=>34,
             'D'=>35,'E'=>36,'F'=>37,'G'=>38,'H'=>39,
             'I'=>40,'J'=>41);

  if(exists($sam2pool{$1})) {

    if($sam2pool{$1} eq "1.8") {
      $headLink = ' ';

      &filterPEreads($in1,$in2,$out1,$out2,$out3,\%asc2ps,$trimQ,$delAQ,$minLen,$Npercent,$headLink);

    }
    elsif($sam2pool{$1} eq "1.3") {
      my (@ins,@outs,@thr);

      ($ins[0],$ins[1]) = ($in1,$in2);
      $in1 =~ s/\.gz$//;
      $in2 =~ s/\.gz$//;
      $outs[0] = $in1."s";       #--it depends--
      $outs[1] = $in2."s";

      for(my $i=0;$i<=1;$i+=1) {
        $thr[$i] = threads->create(\&transQR,$ins[$i],$outs[$i]);
      }

      for(my $i=0;$i<=1;$i+=1) {
        $thr[$i]->join;
      }

      $headLink = '/';
      $in1 = $outs[0].".gz";
      $in2 = $outs[1].".gz";

      &filterPEreads($in1,$in2,$out1,$out2,$out3,\%asc2ps,$trimQ,$delAQ,$minLen,$Npercent,$headLink);

     system("rm $in1 $in2");
    }
  }

}


sub transQR {

  my ($in,$out) = @_;

  #system("perl fq2std.pl sol2std $in $out");

  &fq2sg($in,$out);

}


sub fq2sg {

  my ($in,$out) = @_;

  my @conv_table;
  for (-64..64) {
    $conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
  }

  open(my $FR,"gzip -dc $in |");
  open(my $FW,"| gzip >$out.gz");

  while($_=<$FR>) {
    if (/^@/) {

      print $FW $_;
      $_ = <$FR>; print $FW $_; $_ = <$FR>; $_ = <$FR>;
      chomp;

      my @t = split('', $_);

      my $qual = '';
      $qual .= $conv_table[ord($_)] for (@t);

      print $FW "+\n$qual\n";
    }
  }

  close($FR);
  close($FW);

}


sub filterPEreads {

  my ($in1,$in2,$out1,$out2,$out3,$asc2ps,$trimQ,$delAQ,$minLen,$Npercent,$headLink) = @_;

  open(my $FR1,"gzip -dc $in1 |");
  open(my $FR2,"gzip -dc $in2 |");
  open(my $FW1,"| gzip >$out1.gz");
  open(my $FW2,"| gzip >$out2.gz");
  open(my $FW3,">>$out3");

  my ($numPReadS,$numBaseS,$numPReadQ,$numBaseQ,$numBaseQA) = (0,0,0,0,0,0);
  my (@line1,@line2,$sumQ,$i,$templine,$templineN,$len,$len1,$len2,$Nnum,%read2rec);
  my ($line11,$line12,$line13,$line14,$line21,$line22,$line23,$line24);

  while(defined($line11=<$FR1>)) {

    my ($n1,$n2) = ("a","b");
    while(!($line11 =~ /^@/)) {
      print $line11;
      $line11=<$FR1>;
    }

    if($line11 =~ /^@([^\s]+)($headLink)1/) {   #--it depends------------
      $n1 = $1;
    }
    $line12 = <$FR1>;
    chomp($line12);
    $line13 = <$FR1>;
    $line14 = <$FR1>;
    chomp($line14);

    $line21 = <$FR2>;
    while(!($line21 =~ /^@/)) {
      print $line21;
      $line21=<$FR2>;
    }

    if($line21 =~ /^@([^\s]+)($headLink)2/) {   #--it depends------------
      $n2 = $1;
    }

    $line22 = <$FR2>;
    chomp($line22);
    $line23 = <$FR2>;
    $line24 = <$FR2>;
    chomp($line24);

    #if($n1 ne $n2) {
    #  print "$line11\n$n1\n\n$line21\n$n2\n";
    #  exit;
    #}

    #----remove reads with 5% N and average quality low than required
    my $flag = 1;
    $numPReadS += 1;
    $len1 = length($line12);
    $len2 = length($line22) + $len1;
    $numBaseS += $len2;

    &NlowQ($line12,$line14,$asc2ps,\$flag,$delAQ,$Npercent);

    if($flag==1) {

      &NlowQ($line22,$line24,$asc2ps,\$flag,$delAQ,$Npercent);

      #----check duplicated reads
      if($flag==1) {

        &rmDuplicates($line12,$line22,\%read2rec,\$flag);

        #----trim the 3' bases
        if($flag==1) {
          my ($tlen,$tlenA) = (0,0);
          &trim3tail($line12,$line14,$asc2ps,\$flag,$trimQ,$minLen,\$tlen,\$tlenA);

          if($flag==1) {
            &trim3tail($line22,$line24,$asc2ps,\$flag,$trimQ,$minLen,\$tlen,\$tlenA);

            if($flag==1) {
              $numPReadQ += 1;
              $numBaseQ += $tlen;
              $numBaseQA += $tlenA;

              print $FW1 $line11;
              print $FW1 $line12,"\n";
              print $FW1 $line13;
              print $FW1 $line14,"\n";
              print $FW2 $line21;
              print $FW2 $line22,"\n";
              print $FW2 $line23;
              print $FW2 $line24,"\n";

            }
            else {
              #print "flag 0 in trim read2.\n";
            }
          }
          else {
            #print "flag 0 in trim read1.\n";
          }
        }
        else {
          #print "flag 0 in duplicates.\n";
        }
      }
      else {
        #print "flag 0 in read2 N percent.\n";
      }
    }
    else {
      #print "flag 0 in read1 N percent.\n";
    }
  }

  close($FR1);
  close($FR2);
  close($FW1);
  close($FW2);

  my $ratioRead = $numPReadQ/$numPReadS;
  my $ratioBase = $numBaseQ/$numBaseS;

  print $FW3 "$in1\t$numPReadQ\t$numPReadS\t$ratioRead\t$numBaseQ\t$numBaseS\t$ratioBase\t$numBaseQA\n";

  close($FW3);

}


sub trim3tail {

  my ($seq,$qua,$asc2ps,$flag,$trimQ,$minLen,$tlen,$tlenA) = @_;

  my $line1ta = $seq;
  my $line1tb = $qua;

  my $len = length($line1tb);
  $$tlenA += $len;

  my $temp = chop($line1tb);
  my $tempN = chop($line1ta);

  $len -= 1;

  while($asc2ps->{$temp}<$trimQ && $len>0) {
    $temp = chop($line1tb);
    $tempN = chop($line1ta);
    $len -= 1;
  }

  if($asc2ps->{$temp}>=$trimQ) {
    $line1tb = $line1tb.$temp;
    $line1ta = $line1ta.$tempN;
    $len += 1;
  }

  if($len<$minLen) {
    $$flag = 0;
  }
  else {
    $$tlen += $len;
  }

}


sub rmDuplicates {

  my ($line1,$line2,$read2rec,$flag) = @_;

  my $subtemp1 = $line1."-".$line2;
  my $temp1re = reverse($line1);         #only removing duplicated reads from one cluster
  $temp1re =~ tr/ATCGatcg/TAGCtagc/;

  my $temp2re = reverse($line2);
  $temp2re =~ tr/ATCGatcg/TAGCtagc/;

  my $subtemp2 = $temp2re."-".$temp1re;
  if(!exists($read2rec->{$subtemp1}) && !exists($read2rec->{$subtemp2})){
    $read2rec->{$subtemp1} = "y";
    $read2rec->{$subtemp2} = "y";
  }
  else {
    $$flag = 0;
  }

}


sub NlowQ {

  my ($seq,$qua,$asc2ps,$flag,$delAQ,$Npercent) = @_;

  my $line1ta = $seq;
  my $line1tb = $qua;

  my $len = length($line1tb);
  my $len1 = length($line1tb);

  my $temp = chop($line1tb);
  my $sumQ = $asc2ps->{$temp};

  my $tempN = chop($line1ta);
  my $sumN = 0;

  if($tempN eq "N") {
    $sumN += 1;
  }
  $len -= 1;

  while($len>0) {
    $temp = chop($line1tb);
    $sumQ += $asc2ps->{$temp};
    $tempN = chop($line1ta);

    if($tempN eq "N") {
      $sumN += 1;
    }

    $len -= 1; #length($line1[3]);
  }

  if($sumQ/$len1<$delAQ || $sumN/$len1>$Npercent) {
    $$flag = 0;
  }

}


sub readSam {

  my ($in,$sam2fq) = @_;

  open(my $FR,$in);
  while($_=<$FR>) {
    chomp;

    my @temp = split(/\t/);
    $sam2fq->{$temp[0]} = $temp[1];

  }
  close($FR);

}


sub readFile {

  my ($file,$file2rec,$dir) = @_;

  opendir(DIR,$dir) || die "Error in opening dir $dir\n";
  my $filename;
  while(($filename = readdir(DIR))){

    if($filename=~/^(.+)_\d\.fq.*gz/) {     #it depends--- _hb_1.fq
      $file->{$1} = "Y";
      $file2rec->{$filename} = "Y";
    }
  }

  closedir(DIR);

}



#!/usr/bin/perl -w
use strict;

my $genome = $ARGV[0]; ##- 
my $pe1 = $ARGV[1]; ##- sperated by comma
my $pe2 = $ARGV[2]; ##- sperated by comma

my $cpu = 20;
my ($pwd, $cmdString);
my $HIST2build = '/home/caix/miniconda3/envs/py2.7/bin/hisat2-build';
my $HIST2 = '/home/caix/miniconda3/envs/py2.7/bin/hisat2';
my $SAMTOOLS = '/home/caix/miniconda3/envs/py2.7/bin/samtools';

# Step 2: HISAT2
print STDERR "\n============================================\n";
print STDERR "Step 2: HISAT2\n";
mkdir "2.hisat2" unless -e "2.hisat2";
unless (-e "2.hisat2.ok") {
    chdir "2.hisat2";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    # 构建基因组hisat2索引数据库
    $cmdString = "$HIST2build -p 20  ../$genome genome &> hisat2-build.log\n";
    unless (-e "hisat2-build.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "hisat2-build.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 将RNA-Seq reads比对到参考基因组
    # 能处理单端，双端，链特异性与链非特异性数据
    my $input;
    if ($pe1 && $pe2) {
        $input = "-1 $pe1 -2 $pe2";
    }
    $cmdString = "$HIST2  -x genome -p $cpu $input -S hisat2.sam --min-intronlen 20 --max-intronlen 20000 --dta --score-min L,0.0,-0.4  2> hisat2.log";
    unless (-e "hisat2.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "hisat2.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    $cmdString = "$SAMTOOLS  sort -@ $cpu -o hisat2.sorted.bam -O BAM hisat2.sam";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    $cmdString = "$SAMTOOLS  view -h hisat2.sorted.bam > hisat2.sorted.sam";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    chdir "../";
    open OUT, ">", "2.hisat2.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 2 for the file 2.hisat2.ok exists\n";
}

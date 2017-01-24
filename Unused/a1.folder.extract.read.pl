#!/usr/perl/bin
# run the script: time perl 00.script/a1.folder.extract.read.pl 01.data/13.bowtie/$sample 01.data/14.AS.count/$sample 01.data/00.PriorData/PhQ.transcript.fa

use strict;

my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $fafile = shift @ARGV;
system("mkdir -p $tgtfolder");

my @temp = split(/\//, $srcfolder);
my $sample = pop @temp;

opendir(SRC, $srcfolder);
my @sams = sort(grep(/^SRX/, readdir(SRC)));
closedir(SRC);

system(":> $tgtfolder/$sample.count.txt");
foreach my $sam (@sams){
	system("time perl 00.script/a1.extract.read.number.pl $srcfolder/$sam/bowtie.out.bam $fafile >> $tgtfolder/$sample.count.txt");
}
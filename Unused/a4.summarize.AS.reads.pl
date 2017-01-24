#!/usr/perl/bin
# run the script: time perl 00.script/a4.summarize.AS.reads.pl 01.data/14.AS.count 01.data/00.PriorData 01.data/14.AS.count/as_read_count.txt

use strict;

my $srcfolder = shift @ARGV;
my $reffolder = shift @ARGV;
my $tgtfile = shift @ARGV;

opendir(SRC, $srcfolder);
my @sams = sort(grep(/^[A-Z]\w+/, readdir(SRC)));
closedir SRC;
open(TGT, ">$tgtfile");
print TGT "Experiment\tSample\tSample_Name\tMenA.1\tMenA.2_left\tMenA.2_right\tMenG.1_1\tMenG.1_2\tMenG.1_3\tMenG.2_1_left\tMenG.2_1_right\tMenG.2_2_left\tMenG.2_2_right\tMenG.2_3_left\tMenG.2_3_right\n";

foreach my $sam (@sams){
	print "$sam\n";
	open(REF, "$reffolder/$sam");
	my %hash = ();
	foreach my $line (<REF>){
		chomp $line;
		my @lines = split(/\s+/, $line);
		if(not exists $hash{$lines[1]}){
			$hash{$lines[1]} = $lines[2];
			#print "$lines[1]\t$lines[2]\n";
		}
	}
	close REF;
	
	open(SRC, "$srcfolder/$sam/$sam.count.txt");
	while(my $line = <SRC>){
		my @lines = split(/\s+/, $line);
		my $sample = shift @lines;
		if(exists $hash{$sample}){
			print TGT "$sam\t$sample\t$hash{$sample}\t", join("\t", @lines), "\n";
			#print "$sam\t$sample\t$hash{$sample}\t", join("\t", @lines), "\n";
		}
	}
	close SRC;
}

close TGT;

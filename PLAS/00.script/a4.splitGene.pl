###########################################################
###	@Author		Xi Gu									###
###	@Time		Oct 28, 2014							###
###	@Function	Split Query into Multiple Files			###
###########################################################

#!/usr/bin/perl
# run the script: time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/proteome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.SplitGenes/01.Protein/run.0 1000

use strict;

my $seqFile = shift @ARGV;			# proteome seq file
my $refFile = shift @ARGV;			# gene id file
my $outFolder = shift @ARGV;			# output file
my $splitSize = shift @ARGV;


open(SEQ, $seqFile) or die "Cannot open $seqFile: $!";
						print "Opened seqfile $seqFile\n";
my %hash = ();
my $gene = 0;
						print "Initial gene: $gene\n";
while (my $line = <SEQ>){
	chomp $line;				print "Current line: $line\n";
	$line =~ s/\s+$//g;			print "Post-sed line: $line\n";
	if($line =~ /^>/){			print "I see a gene ID!\n";
		$gene = $line;			print "Gene from line: $gene\n";
		$gene =~ s/>//;			print "Post-sed gene: $gene\n";
		$hash{$gene} = ();		print "Hash for gene: $hash{$gene}\n";
	}else{					print "I do not see a gene ID!\n";
		push @{$hash{$gene}},$line;	print "I do not know what this does: $line\n";
	}
}

open(REF, $refFile) or die "Cannot open $refFile: $!";
my $i = 1;					print "\n\nOpened reference!\n\n";
my $group = $i * $splitSize;		# use splitSize to name the folder and file
system("mkdir -p $outFolder/$group");
open(TGT, ">$outFolder/$group/$group.fasta")or die "Cannot write $group.fasta: $!";

foreach my $line (<REF>){
	chomp $line;				print "Reference line: $line\n";
	my @lines = split(/\t/, $line);
	$line = $lines[0];			print "Post split line: $line\n";
	if($line eq ""){	# come across the blank line, reach the next group
		$i = $i + 1;			print "I see a blank line!\n";
		close TGT;	# close previous file
		$group = $i * $splitSize;
		system("mkdir -p $outFolder/$group");	# make a new file for current group
		open(TGT, ">$outFolder/$group/$group.fasta")or die "Cannot write $group.fasta: $!";
		next;
	}
	
	## not come across with the blank line, so print to the file
	print TGT ">$line\n";
	print TGT join("\n", @{$hash{$line}}), "\n";
}

if(-z "$outFolder/$group/$group.fasta"){
	#system("rm -r $outFolder/$group");
}

close REF;
close SEQ;
close TGT;

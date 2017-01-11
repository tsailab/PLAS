#!/usr/perl/bin
# run the script: time perl 00.script/a1.extract.read.number.pl 13.abundance/PhQ/SRX221648/PhQ.bam 01.data/00.PriorData/PhQ.transcript.fa 01.data/00.PriorData/contig2gene.txt 13.abundance/PhQ/read.count.txt 100

use strict;
use Bio::DB::Sam;

my $bamfile = shift @ARGV;
my $fafile = shift @ARGV;

my @temp = split(/\//, $bamfile);
pop @temp;
my $sample = pop @temp;

# high level API
my $sam = Bio::DB::Sam->new(-bam  =>"$bamfile",
                             -fasta=>"$fafile",
							 -autoindex=>1,
                             );

#my @targets    = $sam->seq_ids;
#open(TGT, ">$tgtfile");

my @genes = ("AT1G60600.1","AT1G60600.2","AT1G23360.1","AT1G23360.2");
my @ranges = ([520],
			  [520, 528],
			  [170, 234, 326],
			  [170, 261, 325, 409, 501, 577]);

## for MenA and MenG
print "$sample\t";
#print "$sample\n";
for my $i(0..$#genes){
	my $gene = $genes[$i];
	my @range = @{$ranges[$i]};
	for my $range (@range){
		my @align = $sam->get_features_by_location(-seq_id => $gene,
                                                -start  => $range,
                                                -end    => $range + 1);
		print scalar(@align), "\t";
		#foreach my $a (@align){
		#	my $seqid  = $a->seq_id;
		#	my $start  = $a->start;
		#	my $end    = $a->end;
		#	my $strand = $a->strand;
		#	print "$seqid\t$start\t$end\t$strand\n";
		#}
	}
}
print "\n";

#close TGT;
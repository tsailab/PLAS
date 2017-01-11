#!/usr/perl/bin
# run the script: time perl 00.script/a1.extract.read.number.pl 13.abundance/PhQ/SRX221648/PhQ.bam 01.data/00.PriorData/PhQ.transcript.fa 01.data/00.PriorData/contig2gene.txt 13.abundance/PhQ/read.count.txt 100

use strict;
use Bio::DB::Sam;

my $bamfile = shift @ARGV;
my $fafile = shift @ARGV;
my $reffile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $winsize = shift @ARGV;

 # high level API
 my $sam = Bio::DB::Sam->new(-bam  =>"$bamfile",
                             -fasta=>"$fafile",
							 -autoindex=>1,
                             );

 #my @targets    = $sam->seq_ids;
 open(REF, $reffile);
 open(TGT, ">$tgtfile");

 <REF>;
 foreach my $line (<REF>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	my $gene = $lines[0];
	my $len = $lines[2];
	my @distn = ();
	print "$gene\n";
	for(my $start=0; $start <= $len-$winsize+1; $start ++){
		my @alignments = $sam->get_features_by_location(-seq_id => $gene,
                                                 -start  => $start,
                                                 -end    => $start + $winsize - 1);
		push @distn, scalar(@alignments);
	}
	print TGT "$lines[0]\t$lines[1]\t", join("\t", @distn), "\n";
	print scalar(@distn), "\n";
	
	#my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$gene);
	#my @data       = $coverage->coverage;
	#print join(",", @data), "\n";
 }
 
 close REF;
 close TGT;
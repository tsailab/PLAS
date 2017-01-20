#!/usr/bin/perl
# time perl 00.script/a5.retrieveseq.GOI.pl 14.annotation/classIII.pero.out 01.data/06.TargetTranscriptome/transcriptome.fa 14.annotation/classIII.pero.seq.fasta
use strict;
use Bio::SeqIO;

my $srcfile = shift @ARGV;
my $seqfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(SRC, $srcfile) or die;
open(TGT, ">$tgtfile") or die;

my %hash =();
foreach my $line (<SRC>){
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[0]}){
		$hash{$lines[0]} = 0;
	}
}

my $seqio = Bio::SeqIO->new(-file => $seqfile, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq){
	my $id = $seq->id;
	if(exists $hash{$id}){
		my $seq = $seq->seq;
		print TGT ">$id\n";
		print TGT "$seq\n";
	}
}

close SRC;
close TGT;
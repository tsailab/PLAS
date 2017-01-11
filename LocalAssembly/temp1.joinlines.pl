#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => 'fasta');
open(TGT, ">$tgtfile");

while(my $seqio = $seqio_obj->next_seq){
    my $id = $seqio->id;
    my $desc = $seqio->description;
    my $seq = $seqio->seq;

    print TGT ">$id $desc\n$seq\n";
}
close TGT;


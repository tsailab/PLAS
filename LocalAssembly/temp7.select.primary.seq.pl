#!/usr/bin/perl
use strict;
use Bio::SeqIO;

my $reffile = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

my %hash = ();
my $refio_obj = Bio::SeqIO->new(-file => $reffile, -format => "fasta");
my $srcio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");
open(TGT, ">$tgtfile");

while(my $seqio = $refio_obj->next_seq){
    my $id = $seqio->id;
    if(not exists $hash{$id}){
        $hash{$id} = 0;
    }
}

while(my $seqio = $srcio_obj->next_seq){
    my $id = $seqio->id;
    if(exists $hash{$id}){
        my $seq = $seqio->seq;
        my @ids = split(/\./, $id);
        pop @ids;
        $id = join(".", @ids);
        print TGT ">$id\n$seq\n";
    }
}

close TGT;

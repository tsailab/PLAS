#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $srcfile = shift @ARGV;
my $seqfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $mode = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

my %hash = ();
while(my $line = <SRC>){
    chomp $line;
    my @lines = split(/\t/, $line);
    if($lines[0] eq $mode){
        if(not exists $hash{$lines[1]}){
            $hash{$lines[1]} = 0;
        }
    }
}
close SRC;

my $seqio_obj = Bio::SeqIO->new(-file => $seqfile, -format => 'fasta');
while(my $seqio = $seqio_obj->next_seq){
    my $id = $seqio->id;
    my $desc = $seqio->description;
    my $seq = $seqio->seq;

    if(exists $hash{$id}){
        print TGT ">$id $desc\n$seq\n";
    }
}
close TGT;

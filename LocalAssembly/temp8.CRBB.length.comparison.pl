#!/usr/bin/perl

use strict;
use Bio::SeqIO;

my $reffile = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $seqio_obj = Bio::SeqIO->new(-file => $reffile, -format => 'fasta');
open(TGT, ">$tgtfile");

my %hash = ();
while(my $seqio = $seqio_obj->next_seq){
    my $id = $seqio->id;
    my $seq = $seqio->seq;
    my $len = length($seq);
    if(not exists $hash{$id}){
        $hash{$id} = $len;
}
}

open(SRC, $srcfile);
<SRC>;
my $q_totallen = 0;
my $hit_totallen = 0;
my $gain = 0;
my $loss = 0;
my %trinity = 0;
my %local = 0;
while(my $line = <SRC>){
    chomp $line;
    my @lines = split(/,/, $line);
    my $hitname = pop @lines;
    if($hitname eq 'NA'){next;}
    my $qname = $lines[0];
    my $qlen = $lines[1];
    my $hitlen = $hash{$hitname};
    
    my $diff = $hitlen - $qlen;
    if(not exists $trinity{$qname}){
        $trinity{$qname} = 0;
        $q_totallen += $qlen;
    }
    if(not exists $local{$hitname}){
        $local{$hitname} = 0;
        $hit_totallen += $hitlen;
    }

    if($diff > 0){
        $gain += $diff;
    }elsif($diff <0){
        $loss += abs($diff);
    }

    print TGT "$qname\t$qlen\t$hitname\t$hitlen\n";
}

print "Query total length: $q_totallen\n";
print "Hit total length: $hit_totallen\n";
print "Total gain length: $gain\n";
print "Total loss length: $loss\n";

close SRC;
close TGT;

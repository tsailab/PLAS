#!/usr/bin/perl -w 
use strict;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

while(my $line = <SRC>){
    chomp $line;
    my @lines = split(/\s+/, $line);
    shift @lines;

    my $i = 0;
    my @genes = ();
    foreach my $w (@lines){
        if($w =~ /^Potri/){
            $i ++;
            push @genes, $w;
        }
    }

    my @out = map {"$_\t$i"} @genes;
    print TGT join("\n", @out), "\n";
}

close SRC;
close TGT;

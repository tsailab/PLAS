#!/usr/bin/perl

use strict;
my $srcfile = shift @ARGV;

my $i = 0;
my $prev = 0;
my $gene = 0;

open(SRC, $srcfile);
while(my $line = <SRC>){
    chomp $line;
    $line =~ s/^>//g;
    my @lines = split(/_/, $line);
    $gene = join("_", @lines[0..1]);
    if($gene ne $prev){
        $i ++;
    }
    $prev = $gene;
} 

print "$i\n";
close SRC;

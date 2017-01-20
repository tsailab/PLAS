#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

## read in parameters required by the script
my $blastfile = shift @ARGV;              ## blast result
my $srcfile = shift @ARGV;               ## assembled contig file
my $tgtfile = shift @ARGV;              ## target folder

## read in seq length
my %query_length = ();

my $query_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");
while(my $seq_obj = $query_obj->next_seq){
    my $id = $seq_obj->id;
    my $seq = $seq_obj->seq;
    my $len = length($seq);
    if(not exists $query_length{$id}){
        $query_length{$id} = $len;
    }
}

open(BLS, $blastfile);
open(TGT, ">$tgtfile");
my %hash = ();

while(my $line = <BLS>){
    chomp $line;
    my @lines = split(/\t/, $line);
    my $qname = $lines[0];
    my $qlen = $query_length{$qname};
    my $hitname = $lines[1];
    my $qstart = $lines[6];
    my $qend = $lines[7];
    my $hitstart = $lines[8];
    my $hitend = $lines[9];
    my $evalue = $lines[10];
    my $bits = $lines[11];
    my $frame = 0;

    ## find open reading frame
    if($qstart < $qend and $hitstart < $hitend){
        if($qstart % 3 == 1){
            $frame = 1;
        }elsif($qstart % 3 == 2){
            $frame = 2;
        }else{
            $frame = 3;
        }
    }elsif($qstart > $qend and $hitstart < $hitend){
        if(($qlen - $qstart) % 3 == 0){
            $frame = -1;
        }elsif(($qlen - $qstart) % 3 == 1){
            $frame = -2;
        }else{
            $frame = -3;
        }
    }

    if(not exists $hash{$qname}){
        $hash{$qname} = $frame;
    }    
}

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");
while(my $seq_obj = $seqio_obj->next_seq){
    my $id = $seq_obj->id;
    my $desc = $seq_obj->description;
    my $seq = $seq_obj->seq;
    
    if(exists $hash{$id}){
        my $proseq = &Find_ORF($seq_obj, $hash{$id});
        print TGT ">$id\n";
        print TGT "$proseq\n";
    }
}

close TGT;

########################
sub Find_ORF{
    my $seq_obj = shift @_;
    my $frame = shift @_;
    my $pro_seq = 0;
    my $protein = 0;

    ## get all six posssible translation
    if($frame == 1){$pro_seq = $seq_obj->translate(-frame => 0)->seq;}
    elsif($frame == 2){$pro_seq = $seq_obj->translate(-frame => 1)->seq;}
    elsif($frame == 3){$pro_seq = $seq_obj->translate(-frame => 2)->seq;}
    elsif($frame == -1){$pro_seq = $seq_obj->revcom->translate(-frame => 0)->seq;}
    elsif($frame == -2){$pro_seq = $seq_obj->revcom->translate(-frame => 1)->seq;}
    elsif($frame == -3){$pro_seq = $seq_obj->revcom->translate(-frame => 2)->seq;}

    ## return value
    return $pro_seq;
}


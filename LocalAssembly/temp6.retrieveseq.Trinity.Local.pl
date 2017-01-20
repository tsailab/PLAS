#!/usr/bin/perl
# run the script: time perl 00.script/temp6.retrieveseq.Trinity.Local.pl 12.Local.Trinity.comparison/Trinity_Local_intersect.txt 13.Local.Trinity.alignment 01.data/06.TargetTranscriptome/transcriptome.v1.fa 05.Trinity/Trinity.new.fasta 01.data/00.PriorData/ptr.transcriptome.fa

use strict;
use Bio::SeqIO;

my $srcfile = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $local_nucl = shift @ARGV;
my $trinity_nucl = shift @ARGV;
my $ref_nucl = shift @ARGV;

my %local_seq = ();
my %trinity_seq = ();
my %ref_seq = ();

system("mkdir -p $tgtfolder/01.data");

my $local_nucl_obj = Bio::SeqIO->new(-file => $local_nucl, -format => 'fasta');
my $trinity_nucl_obj = Bio::SeqIO->new(-file => $trinity_nucl, -format => 'fasta');
my $ref_nucl_obj = Bio::SeqIO->new(-file => $ref_nucl, -format => 'fasta');

while(my $seqio = $local_nucl_obj->next_seq){
	my $id = $seqio->id;
	if(not exists $local_seq{$id}{'nucl'}){
		$local_seq{$id}{'nucl'} = $seqio;
	}
}

while(my $seqio = $trinity_nucl_obj->next_seq){
	my $id = $seqio->id;
	if(not exists $trinity_seq{$id}{'nucl'}){
		$trinity_seq{$id}{'nucl'} = $seqio;
	}
}

while(my $seqio = $ref_nucl_obj->next_seq){
	my $id = $seqio->id;
	if(not exists $ref_seq{$id}{'nucl'}){
		$ref_seq{$id}{'nucl'} = $seqio;
	}
}

my $i = 0;
open(SRC, $srcfile);
while(my $line = <SRC>){
	chomp $line;
	if($line =~ /^>/){
        $i ++;
		$line =~ s/>//g;
        my @lines = split(/\s+/, $line);
        $line = $lines[0];
		if($i > 1){close NUC; }
		open(NUC, ">$tgtfolder/01.data/$i.$line.nucl.fasta");
        
        if(exists $ref_seq{$line}){
            print NUC ">$line\n", $ref_seq{$line}{'nucl'}->seq, "\n";
        }else{
            print "Warning: Not existing $line in reference\n";
        }
	}else{
		my @lines = split(/\s+/, $line);
        my $start = $lines[7];
        my $end = $lines[8];
        my $seq = 0;
        if($lines[0] eq 'Trinity'){
		    if(exists $trinity_seq{$lines[1]}){
                if($start < $end){
                    $seq = $trinity_seq{$lines[1]}{'nucl'}->seq;
                }else{
                    $seq = $trinity_seq{$lines[1]}{'nucl'}->revcom->seq;
                }
			    print NUC ">Trinity.$lines[1]\n$seq\n";
		    }else{
			    print "Warning: Not existing $lines[1] in Trinity\n";
		    }
        }elsif($lines[0] eq 'Local'){
		    if(exists $local_seq{$lines[1]}){
                if($start < $end){
                    $seq = $local_seq{$lines[1]}{'nucl'}->seq;
                }else{
                    $seq = $local_seq{$lines[1]}{'nucl'}->revcom->seq;
                }
			    print NUC ">Local.$lines[1]\n$seq\n";
		    }else{
			    print "Warning: Not existing $lines[1] in Local\n";
		    }
        }
	}
}

close SRC;
close NUC;


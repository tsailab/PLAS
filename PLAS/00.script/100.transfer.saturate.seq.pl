#!/usr/bin/perl -w
# run the script: time perl 00.script/10.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.0 07.map.back/03.bowtie.nucl/run.0 01.data/05.SplitGenes/02.Transcript/run.1 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 

use strict;
use Bio::SeqIO;

## read in parameters required by the script
my $srcfolder = shift @ARGV;			## assembled contig folder
my $blastfolder = shift @ARGV;			## blast result folder
my $tgtfolder = shift @ARGV;			## target folder to put contigs for new run
my $reffile = shift @ARGV;				## reference file with gene ID and protein length info
my $mode = shift @ARGV;					## abs: absolute value; pct: percent value
my $cutoff = shift @ARGV;				## absolute AA number, or percent
my $sleeptime = shift @ARGV;
my ($run) = $srcfolder =~ /\/run\.([0-9]+)/;

## start running the script

system("mkdir -p $tgtfolder");

opendir(SRC, $blastfolder) or die "ERROR: Cannot open $blastfolder: $!";
my @subs = sort(grep(/^[0-9]+$/, readdir SRC));
closedir SRC;

foreach my $sub (@subs){
	
	my $trinity_obj = Bio::SeqIO->new(-file => "$srcfolder/$sub.trinity/Trinity.new.fasta", -format => 'fasta');
	my %query_length = ();
	while(my $seqio = $trinity_obj->next_seq){
		my $id = $seqio->id;
		my $seq = $seqio->seq;
		my $len = length($seq);
		if(not exists $query_length{$id}){
			$query_length{$id} = $len;
		}
	}
	
	open(BLS, "$blastfolder/$sub/$sub.contigs.blast.out") or die "ERROR: Cannot open $blastfolder/$sub/$sub.contigs.blast.out: $!";
	my %valid = ();
	foreach my $line (<BLS>){
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
		
		if(not exists $valid{$qname} or ($bits > $valid{$qname}->{bits})){
            $valid{$qname} = {bits => $bits, frame => $frame};
        }
		
	}
	close(BLS);
	
	$trinity_obj = Bio::SeqIO->new(-file => "$srcfolder/$sub.trinity/Trinity.new.fasta", -format => 'fasta');
	system("mkdir -p $tgtfolder/$sub");
	open(TGT1, ">$tgtfolder/$sub/$sub.fasta") or die "ERROR: Cannot write $tgtfolder/$sub/$sub.fasta: $!";
	open(TGT2, ">$tgtfolder/$sub/$sub.prot.fasta") or die "ERROR: Cannot write $tgtfolder/$sub/$sub.prot.fasta: $!";
	
	while(my $seqio = $trinity_obj->next_seq){
		my $id = $seqio->id;
		my @desc = split(/\s+/, $seqio->description);
		my $seq = $seqio->seq;
		
		if(exists $valid{$id}){
			print TGT1 ">$id ", $desc[0], "\n";
			print TGT1 "$seq\n";
			
			my ($fullseq, $proseq, $qlen) = &Find_ORF($seqio, $valid{$id}->{frame});
			print TGT2 ">$id $qlen\n";
			print TGT2 "$fullseq\n";
		}
	}
}

close TGT1;
close TGT2;

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

	## find the longest ORF between two stop codons
	my $prolen = 0;
	my $stop = 0;
	my @orfs = split(/\*/, $pro_seq);
	if(scalar @orfs > 1){$stop = 1;}
	foreach my $orf (@orfs){
		if(length($orf) > $prolen){
			$protein = $orf;
			$prolen = length($orf);
		}
	}

	## chop AAs before the start codon
	$protein =~ s/^[A-L|N-Z]+M/M/;
	## add stop codon to the end
	if($stop){$protein = $protein."*";}
	## update protein length
	$prolen = length($protein);
	
	## return value
	return ($pro_seq, $protein, $prolen);
}


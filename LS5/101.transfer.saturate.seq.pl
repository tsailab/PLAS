#!/usr/bin/perl -w
# run the script: time perl 00.script/101.transfer.saturate.seq.pl 10.unmapped.reads.trinity/Trinity.new.fasta 10.unmapped.reads.trinity/blastx.out 10.unmapped.reads.trinity 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 

use strict;

use Bio::SeqIO;
use Bio::SearchIO;

## read in parameters required by the script
my $usage = "usage: $0 blast+.outfmt6 query.fasta db.fasta output_folder [output_prefix=NameOfBlastFileHere] [mode=abs|pct] [cutoff=NumericCutoff]\n\n";

my $blastfile = shift @ARGV or die $usage;              ## blast result
my $srcfile = shift @ARGV or die $usage;				## assembled contig file
my $dbfile = shift @ARGV or die $usage;
my $tgtfolder = shift @ARGV or die $usage;              ## target folder
my $output_prefix = shift @ARGV || "$blastfile";
my $mode = shift @ARGV || "pct";					## abs: absolute value; pct: percent value
my $cutoff = shift @ARGV || "0.02";                 ## absolute AA number, or percent

system("mkdir -p $tgtfolder");
my $logfile = "$tgtfolder/detect.full.length.log";

open(TGT1, ">$tgtfolder/$output_prefix.full.contigs.nucl.fasta") or die "ERROR: Cannot write $tgtfolder/$output_prefix.full.length.contigs.fasta: $!";
open(TGT2, ">$tgtfolder/$output_prefix.incomplete.contigs.nucl.fasta") or die "ERROR: Cannot write $tgtfolder/$output_prefix.incomplete.fasta: $!";
open(TGT3, ">$tgtfolder/$output_prefix.unmapped.contigs.nucl.fasta");
open(PRO1, ">$tgtfolder/$output_prefix.full.contigs.prot.fasta") or die "ERROR: Cannot write $tgtfolder/$output_prefix.full.length.contigs.fasta: $!";
open(PRO2, ">$tgtfolder/$output_prefix.incomplete.contigs.prot.fasta") or die "ERROR: Cannot write $tgtfolder/$output_prefix.incomplete.fasta: $!";

open(BLS, $blastfile);

## read in seq length
my %query_length = ();
my %db_length = ();

my $query_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");
while(my $seq_obj = $query_obj->next_seq){
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	my $len = length($seq);
	if(not exists $query_length{$id}){
		$query_length{$id} = $len;
	}
}
my $db_obj = Bio::SeqIO->new(-file => $dbfile, -format => "fasta");
while(my $seq_obj = $db_obj->next_seq){
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	my $len = length($seq);
	if(not exists $db_length{$id}){
		$db_length{$id} = $len;
	}
}

## read in blast results
my %full_length = ();
my %incomplete_length = ();

while(my $line = <BLS>){
	chomp $line;
	my @lines = split(/\t/, $line);
	my $qname = $lines[0];
	my $qlen = $query_length{$qname};
	my $hitname = $lines[1];
	my $hitlen = $db_length{$hitname};
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
	
	my $left = ($hitstart, $hitend)[$hitstart > $hitend];		## subject start
	my $right = ($hitstart, $hitend)[$hitstart < $hitend];		## subject end

	if($left == 1 and $right == $hitlen){
		if(not exists $full_length{$qname} or ($bits > $full_length{$qname}->{bits})){
			$full_length{$qname} = { query_id => $qname,
									 hit_id => $hitname,
									 query_len => $qlen,
									 hit_len => $hitlen,
									 query_start => $qstart,
									 query_end => $qend,
									 hit_start => $hitstart,
									 hit_end => $hitend,
									 evalue => $evalue,
									 bits => $bits,
									 frame => $frame
			};
			#print "$qname\t$qlen\t$hitname\t$hitlen\t$qstart\t$qend\t$frame\n";
			if(exists $full_length{$qname}){
				#print "Warnings: $qname is mapped to the same gene with multiple hit.\n";
			}
		}else{
			#print "Warnings: $qname is mapped to the same gene with multiple hit.\n";
		}
		
	}else{
		if(not exists $incomplete_length{$qname} or ($bits > $incomplete_length{$qname}->{bits})){
			$incomplete_length{$qname} = { query_id => $qname,
										   hit_id => $hitname,
										   query_len => $qlen,
										   hit_len => $hitlen,
										   query_start => $qstart,
										   query_end => $qend,
										   hit_start => $hitstart,
										   hit_end => $hitend,
										   evalue => $evalue,
										   bits => $bits,
										   frame => $frame
										   };
			#print "$qname\t$qlen\t$hitname\t$hitlen\t$qstart\t$qend\t$frame\n";
			if(exists $incomplete_length{$qname}){
				#print "Warnings: $qname is mapped to the same gene with multiple hit.\n";
			}
		}else{
			#print "Warnings: $qname is mapped to the same gene with multiple hit.\n";
		}
	}        
}

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");

while(my $seq_obj = $seqio_obj->next_seq){
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	my @temp = split(/\s+/, $seq_obj->desc);
	my $seqlen = $temp[0];
	
	if(exists $full_length{$id}){
		my ($proseq, $qlen) = &Find_ORF($seq_obj, $full_length{$id}->{frame});
		print TGT1 ">$id len=", $full_length{$id}->{query_len}, " ", 
					$full_length{$id}->{hit_id}, " ",
					$full_length{$id}->{hit_len}, " ", 
					$qlen, " blast_full ", 
					$full_length{$id}->{frame}, "\n";
		print TGT1 "$seq\n";
		print PRO1 ">$id len=", $full_length{$id}->{query_len}, " ", 
					$full_length{$id}->{hit_id}, " ",
					$full_length{$id}->{hit_len}, " ", 
					$qlen, " blast_full ", 
					$full_length{$id}->{frame}, "\n";
		print PRO1 "$proseq\n";
	}elsif(exists $incomplete_length{$id}){
		my $hitlen = $incomplete_length{$id}->{hit_len};
		my ($proseq, $qlen) = &Find_ORF($seq_obj, $incomplete_length{$id}->{frame});
		
		if(($mode eq "abs" and $qlen > ($hitlen-$cutoff)) || ($mode eq "pct" and ($qlen/$hitlen) > (1-$cutoff))){
				## consider as full assembled contig
				print TGT1 ">$id len=", $incomplete_length{$id}->{query_len}, " ", 
						$incomplete_length{$id}->{hit_id}, " ",
						$incomplete_length{$id}->{hit_len}, " ", 
						$qlen, " infer_full ", 
						$incomplete_length{$id}->{frame}, "\n";
				print TGT1 "$seq\n";
				print PRO1 ">$id len=", $incomplete_length{$id}->{query_len}, " ", 
						$incomplete_length{$id}->{hit_id}, " ",
						$incomplete_length{$id}->{hit_len}, " ", 
						$qlen, " infer_full ", 
						$incomplete_length{$id}->{frame}, "\n";
				print PRO1 "$proseq\n";
		}elsif($mode eq "abs" or $mode eq "pct"){
				print TGT2 ">$id len=", $incomplete_length{$id}->{query_len}, " ", 
						$incomplete_length{$id}->{hit_id}, " ",
						$incomplete_length{$id}->{hit_len}, " ", 
						$qlen, " ", 
						$incomplete_length{$id}->{frame}, "\n";
				print TGT2 "$seq\n";
				print PRO2 ">$id len=", $incomplete_length{$id}->{query_len}, " ", 
						$incomplete_length{$id}->{hit_id}, " ",
						$incomplete_length{$id}->{hit_len}, " ", 
						$qlen, " ", 
						$incomplete_length{$id}->{frame}, "\n";
				print PRO2 "$proseq\n";
		}else{
			print "ERROR: Please specify correct mode: abs or pct\n";
			die;
		} ## mode if else
	}else{
		print TGT3 ">$id $seqlen\n";
		print TGT3 "$seq\n";
	}
}

close TGT1;
close TGT2;
close TGT3;
close PRO1;
close PRO2;

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
	return ($pro_seq, $prolen);
}




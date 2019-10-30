#!/usr/bin/perl -w
# run the script: time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/02.Transcript/01.blast/ncbi.blast.all.out 01.data/03.MCL/02.Transcript/02.mcl/mcl.graph.txt ncbi

use strict;

my $srcfile = shift @ARGV;	# blast output as input
my $tgtfile = shift @ARGV;	# outputfile

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

my $preq = "";
my $pres = "";
my $score = 0.0;
my $count = 0;
my $average = 0.0;
my $iter = 0;
my $e = 0;

foreach my $line(<SRC>){
	$iter ++;
	chomp $line;
	my @lines = split(/\t/, $line);
	my ($q, $s) = @lines[0..1];
	$e = $lines[10];
	
	## take log of evalue
	if($e == 0.0){
		$e = 181;
	}else{
		$e = -&log10($e);
	}
	
	## get gene id from transcript id
	my @temp = split(/\./, $q);
	#pop @temp;
	$q = join(".", @temp);
	@temp = split(/\./, $s);
	#pop @temp;
	$s = join(".", @temp);
	#print "$q\t$s\n";
	
	## combine multiple segment between the same gene pair
	if($q eq $preq and $s eq $pres){
		$score += $e;
		$count += 1;
	}else{
		if($iter == 1){
			$score = $e;
			$count = 1;
			$preq = $q;
			$pres = $s;
			next;
		}
		if($count == 0){
			$average = $score;
		}else{
			$average = $score/$count;	## take average of multiple segments
		}
		
		## print to the output file
		print TGT  "$preq\t$pres\t$average\n";
		$score = $e;
		$count = 1;
		$preq = $q;
		$pres = $s;
	}
}

if($count == 0){
	$average = $score;
}else{
	$average = $score/$count;
}
print TGT  "$preq\t$pres\t$average\n";

close(SRC);
close(TGT);

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}


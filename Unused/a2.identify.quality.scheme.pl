#!/usr/perl/bin
# run the script: time perl 00.script/a2.identify.quality.scheme.pl 01.data/00.PriorData/SRA174408.v0 01.data/01.Fastq/SRA174408 01.data/00.PriorData/SRA174408.v1 2

use strict;
system("echo 'Running a2.identify.quality.scheme.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffile = shift @ARGV;
my $srcfolder = shift @ARGV;
my $tgtfile = shift @ARGV;
my $col = shift @ARGV;
my $sleeptime = shift @ARGV;
$col = $col - 1;

my @temp = split(/\//, $reffile);
my $sample = pop @temp;
my $dir = "00.script/00.download.script/$sample";

open(TGT ,">$tgtfile");

## check if previous step has succesfully finished
open(CHK, $reffile);
my @chks = ();
foreach my $line(<CHK>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	my $sam = $lines[$col];
	push @chks, $sam;
}
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	my @stderr = glob("download.$sample.sh.o*");
	my @stdout = glob("download.$sample.sh.e*");
	my @log = glob("$dir/done.*");
	if(!($stderr[0] and $stdout[0] and $log[0])){
		last; # job hasn't finished, no log file
	}else{
		system("grep -E 'ERROR|Error|error' download.$sample.sh.e* >> $dir/summary.error.log\n");
	}
		if(!-s "$dir/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/$chk.fastq") and (!(-s "$srcfolder/$chk/$chk.R1.fastq") or !(-s "$srcfolder/$chk/$chk.R2.fastq"))){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
				}
				else{
					$i++;
				}
			}
			if($i == scalar @chks){
				system("echo 'All jobs have been finished successfully' >> job.monitor.txt"); # error file is empty
				last;
			}
		}else{
			die "ERROR: something went wrong in previous steps\n";	
		}
	sleep $sleeptime;
}
closedir CHK;

open(SRC, $reffile);
my $mode = "";
my @temp = @chks;
foreach my $line(<SRC>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	my $sam = $lines[$col];
	opendir(SUB, "$srcfolder/$sam");
	my @subs = sort(grep(/\w+/, readdir(SUB)));
	closedir SUB;
	
	my $mark = 0;
	my $sub = shift @subs;
	print "Processing $srcfolder/$sam/$sub\n";
	open(SUB, "$srcfolder/$sam/$sub");
	
	my $count = 0;
	while(my $subline = <SUB>){
		$count ++;
		<SUB>; 
		$subline = <SUB>; 
		if(!($subline =~ /^\+/)){
			die "Error: The fastq file is not well-formated.\n";
		}
		$subline = <SUB>;
		chomp $subline;
		if($subline =~ /[!"#\$%&'()*+,-.\/0123456789:;<=>?BCDEFGHIJKLMNOPQRSTUVWXYZ]/){
			if($subline =~ /k/){
				$mode = 'Illumina1.8';
				last;
			}
			$mode = 'Sanger';
		}else{
			if($subline =~ /\\|\]|\^|_|`/){
				$mode = 'Solexa';
				last;
			}
			$mode = 'Illumina1.3';
		}
		
		if($count == 2000){last;}
	}
	
	close SUB;
	print "The format is $mode\n";
	print TGT "$line\t$mode\n";
}

close SRC;
close TGT;

system("echo 'Finished a2.identify.quality.scheme.pl!' >> job.monitor.txt");

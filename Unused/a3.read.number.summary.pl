#!/usr/bin/perl -w
# run the script: time perl 00.script/a3.read.number.summary.pl 01.data/01.Fastq/SRA174408 01.data/03.CutAdapter 01.data/06.ncRNA 01.data/07.cleanRNA Sapelo single-end 10

use strict;
system("echo 'Running a3.read.number.summary.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $cutfolder = shift @ARGV;
my $ncfolder = shift @ARGV;
my $trimfolder = shift @ARGV;
my $mode = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;

my @temp = split(/\//, $reffolder);
my $sample = pop @temp;
my $prev = "00.script/05.trimmomatic.script/$sample";

## check if previous step has succesfully finished
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^\w+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("trimmomatic.$sample.$chk.sh.o*");
		my @stdout = glob("trimmomatic.$sample.$chk.sh.e*");
		my @log = glob("$prev/$chk.done.*");
		if(!($stderr[0] and $stdout[0] and $log[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' trimmomatic.$sample.$chk.sh.e* >> 00.script/05.trimmomatic.script/summary.error.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/05.trimmomatic.script/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$trimfolder/$chk/$chk.fastq") and (!(-s "$trimfolder/$chk/$chk.R1.fastq") or !(-s "$trimfolder/$chk/$chk.R2.fastq"))){
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
	}
	sleep $sleeptime;
}
close CHK;

## start running the script
opendir(SRC,"$reffolder") || die"Cannot open file:$!";
my @sams = sort(grep(/\w+/, readdir(SRC)));
closedir SRC;

open(TGT, ">$trimfolder/read_count_tab");
foreach my $sam(@sams)
{
	print "Processing $sample/$sam\n";
	if($mode eq 'single-end'){
		my $test = "$reffolder/$sam/$sam.fastq";
		my $c1 = &count($test);
		my $c2 = &count("$cutfolder/$sam/$sam.fastq");
		my $c3 = &count("$ncfolder/$sam/$sam.fastq");
		my $c4 = &count("$trimfolder/$sam/$sam.fastq");
		print TGT "$sample\t$sam\t$c1\t$c2\t$c3\t$c4\n";
		print "$sample\t$sam\t$c1\t$c2\t$c3\t$c4\n";
	}elsif($mode eq 'paired-end'){
		my @R = ('R1', 'R2');
		foreach my $R (@R){
			my $c1 = &count("$reffolder/$sam/$sam.$R.fastq");
			my $c2 = &count("$cutfolder/$sam/$sam.$R.fastq");
			my $c3 = &count("$ncfolder/$sam/$sam.$R.fastq");
			my $c4 = &count("$trimfolder/$sam/$sam.$R.fastq");
			print TGT "$sample\t$sam.$R\t$c1\t$c2\t$c3\t$c4\n";
			print "$sample\t$sam.$R\t$c1\t$c2\t$c3\t$c4\n";
		}
	}
}
close TGT;

system("echo 'Finished a3.read.number.summary.pl!' >> job.monitor.txt");

## read count function for fastq
sub count{
	my $file = shift @_;
	open(SAM, $file);
	my $c = 0;
	while(my $line = <SAM>){
		$c ++;
		<SAM>; <SAM>; <SAM>;
	}
	return $c;
}

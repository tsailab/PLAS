#!/usr/bin/perl -w
#running command: time perl 00.script/13.folder.bowtie.pl 01.data/01.Fastq/SRA174408 01.data/13.bowtie/SRA174408 Sapelo single-end 10

use strict;
system("echo 'Running 06.bowtie.folder.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $srcfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $mode = shift @ARGV;
my $thread = 1;

my $currstep = "deseq";
my $prevstep = "retrieve";
my @temp = split(/\//, $reffolder);
my $sample = pop @temp;
my $dir = "00.script/14.$currstep.script/$sample";

## check if previous step has succesfully finished
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^\w+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("$prevstep.$chk.o*");
		my @stdout = glob("$prevstep.$chk.e*");
		if(!($stderr[0] and $stdout[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' $prevstep.$chk.e* >> 00.script/13.$prevstep.script/summary.error.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/13.$prevstep.script/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/bowtie.out.sam")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					#system("rm -f $prevstep.$chk.*");
					system("echo 'Resubmitting the job: $prevstep.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/13.$prevstep.script/$prevstep.$chk.sh");
						#last; # There are some jobs failed
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
opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @sams = sort(grep(/^\w.+/, readdir(SRC)));
closedir SRC;

system("mv $prevstep.* 00.script/13.$prevstep.script/$sample/");
system("rm -rf $dir");
system("mkdir -p $dir");

foreach my $sam (@sams){
	my $shell = "$dir/$currstep.$sam.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N $currstep.$sam\n";
		print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=20gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
		print SHL "module load samtools/1.2\n\n";
		
		print SHL "time samtools view -h grep $srcfolder/$sam/bowtie.bam | grep \"\" > $srcfolder/$sam/PhQ.sam\n";
		print SHL "time samtools view -bS $srcfolder/$sam/PhQ.sam > $srcfolder/$sam/PhQ.bam\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n\n";

		print SHL "time /usr/local/samtools/1.2/bin/samtools view -h grep $srcfolder/$sam/bowtie.bam | grep \"\" > $srcfolder/$sam/PhQ.sam\n";
		print SHL "time /usr/local/samtools/1.2/bin/samtools view -bS $srcfolder/$sam/PhQ.sam > $srcfolder/$sam/PhQ.bam\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	
	close(SHL);
	system("chmod u+x $shell");
	if($platform eq "sapelo"){
    	system("qsub $shell");
	}elsif($platform eq "zcluster"){
		system("qsub -q rcc-30d -pe thread $thread $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
}

system("echo 'Finished 06.bowtie.folder.pl!' >> job.monitor.txt");


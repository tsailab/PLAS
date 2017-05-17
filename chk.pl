## check if previous step has succesfully finished
my $reffolder = "01.data/05.SplitGenes/01.Protein/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("makebowtiedb.$chk.sh.o*");
		my @stdout = glob("makebowtiedb.$chk.sh.e*");
		my @log = glob("00.script/02.makebowtiedb.script/run.$run/$chk.done.*");
		if(!($stderr[0] and $stdout[0] and $log[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' makebowtiedb.$chk.sh.e* >> 00.script/02.makebowtiedb.script/run.$run/summary.error.log\n");
			#system("echo 'success' > 00.script/shell.script/makebowtiedb.$chk.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/02.makebowtiedb.script/run.$run/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$dbfolder/$chk/$chk.1.bt2" and -s "$dbfolder/$chk/$chk.rev.1.bt2")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					#system("rm -f makebowtiedb.$chk.*");
					#system("rm -f 00.script/02.makebowtiedb.script/run.$run/$chk.done.log");
					system("echo 'Resubmitting the job: makebowtiedb.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/02.makebowtiedb.script/run.$run/makebowtiedb.$chk.sh");
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


while(1){
	my $count = 0;
	my @temp = @chks;
	my @temp2 = @sams;
	my $i = 0;
	while(my $chk = shift @temp){
		@temp2 = @sams;
		while(my $sam = shift @temp2){
			my @stderr = glob("bowtie.$seqtype.$chk.$sam.sh.o*");
			my @stdout = glob("bowtie.$seqtype.$chk.$sam.sh.e*");
			my @log = glob("00.script/$logfolder/$chk.$sam.done.*");
			if(!($stderr[0] and $stdout[0] and $log[0])){
				last; # job hasn't finished, no log file
			}else{
				system("grep -E 'ERROR|Error|error' bowtie.$seqtype.$chk.$sam.sh.e* >> 00.script/$logfolder/summary.error.log\n");
				$count ++; # job has finished, add one count
			}
		}			
	}
	@temp = @chks;
	if($count == (scalar @chks) * (scalar @sams)){ # all jobs have finished
		if(!-s "00.script/$logfolder/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
				system("echo 'All jobs have been finished successfully' >> job.monitor.txt"); # error file is empty
				last;
		}else{
			#system("rm -f 00.script/$logfolder/*done.log");
			die "ERROR: something went wrong in previous steps\n";	
		}
	}
	sleep $sleeptime;
}

close CHK;
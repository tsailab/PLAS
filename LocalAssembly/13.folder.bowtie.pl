#!/usr/bin/perl -w
#running command: time perl 00.script/13.folder.bowtie.pl 01.data/01.Fastq/SRA174408 01.data/07.cleanRNA 01.data/00.PriorData/transcriptome.isoform.fa 01.data/13.bowtie/SRA174408 Sapelo single-end

use strict;
system("echo 'Running 06.bowtie.folder.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $qryfolder = shift @ARGV;
my $dbfile = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $mode = shift @ARGV;
my $thread = 2;
my $currstep = "bowtie";
my $dir = "00.script/13.$currstep.script";

## start running the script
opendir(QRY, $qryfolder) or die "ERROR: Cannot open $qryfolder: $!";
my @sams = sort(grep(/^S\w.+/, readdir(QRY)));
closedir QRY;

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
		print SHL "#PBS -l mem=20gb\n";
		print SHL "\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	my $t = $thread / 2;
	if($platform eq "sapelo"){
    	print SHL "module load bowtie2/2.2.4\n";
		print SHL "module load samtools/1.2\n";
	}elsif($platform eq "zcluster"){
		print SHL "PATH=/usr/local/bowtie2/2.2.3/bin/:/usr/local/samtools/1.2/bin/:\$PATH\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	print SHL "mkdir -p $tgtfolder/$sam\n";
			
	if($mode eq "single-end"){
		print SHL "time bowtie2 -q -x $dbfile --no-unal -p $thread -U $qryfolder/$sam/$sam.fastq -S $tgtfolder/$sam/bowtie.out.sam\n";
	}elsif($mode eq "paired-end"){
		print SHL "time bowtie2 -q -x $dbfile --no-unal -p $thread -1 $qryfolder/$sam/$sam.R1.fastq -2 $qryfolder/$sam/$sam.R2.fastq -U $qryfolder/$sam/$sam.R1.U.fastq,$qryfolder/$sam/$sam.R2.U.fastq -S $tgtfolder/$sam/bowtie.out.sam\n";
	}else{
		die "Error: please specify mode as either 'single-end' or 'paired-end'";
	}
	print SHL "samtools view -b -o $tgtfolder/$sam/bowtie.out.bam $tgtfolder/$sam/bowtie.out.sam\n";
	print SHL "rm $tgtfolder/$sam/bowtie.out.sam\n";
	
	close(SHL);
	system("chmod u+x $shell");
	if($platform eq "sapelo"){
    	system("qsub $shell");
	}elsif($platform eq "zcluster"){
		system("qsub -q rcc-30d -pe thread $t $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
}

system("echo 'Finished 06.bowtie.folder.pl!' >> job.monitor.txt");


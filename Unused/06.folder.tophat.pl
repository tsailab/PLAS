#!/usr/bin/perl -w
#running command: time perl 00.script/06.folder.tophat.pl 01.data/01.Fastq/SRA174408 01.data/07.cleanRNA 01.data/00.PriorData/Athaliana_167_TAIR10 01.data/09.tophat single-end Sapelo 10

use strict;
system("echo 'Running 06.folder.tophat.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $srcfolder = shift @ARGV;
my $genome = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $mode = lc(shift @ARGV);
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my $special = shift @ARGV;
my $thread = 4;

my $dir = "00.script/06.tophat.script";

## start running the script
opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @sams = sort(grep(/\w.+/, readdir(SRC)));
closedir SRC;

system("rm -rf $dir");
system("mkdir -p $dir");

## get seq quality scheme
open(SPL, $special);
my %quality = ();
foreach my $line (<SPL>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $quality{$lines[1]}){
		$quality{$lines[1]} = lc($lines[3]);
	}
}
close SPL;

## construct index for reference
=pod
if($platform eq 'sapelo'){
	system("time /usr/local/apps/bowtie2/2.2.4/bin/bowtie2-build -f -q $genome.fa $genome");
}elsif($platform eq 'zcluster'){
	system("time /usr/local/bowtie2/2.2.3/bin/bowtie2-build -f -q $genome.fa $genome");
}
=cut

foreach my $sam (@sams){
	my $shell = "$dir/tophat.$sam.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N tophat.$sam.sh\n";
		print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=20gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
    	print SHL "module load tophat/2.0.13\n\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n\n";
		print SHL "export PATH=/usr/local/samtools/0.1.19/:\${PATH}\n";
		print SHL "export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:\${LD_LIBRARY_PATH}\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	my $t = $thread / 2;
	my $q = '--solexa1.3-quals';
	if(exists $quality{$sam}){
		if($quality{$sam} eq 'sanger' or $quality{$sam} eq 'illumina1.8'){
			$q = "--solexa-quals";
		}
	}else{
		die "Error: quality scheme has not been identified for $sam\n";
	}

	my @R = ();
	if($mode eq "single-end"){
		@R = ("$srcfolder/$sam/$sam.fastq");
	}elsif($mode eq "paired-end"){
		@R = ("$srcfolder/$sam/$sam.R1.fastq_pairs_R1.fastq", "$srcfolder/$sam/$sam.R2.fastq_pairs_R2.fastq");
	}else{
		die "Error: please specify mode as either 'single-end' or 'paired-end'";
	}
	
	if($platform eq "sapelo"){
		print SHL "mkdir -p $tgtfolder/$sam\n";
		print SHL "time tophat -i 30 -I 10000 -g 1 -F 0.05 $q --min-segment-intron 30 --max-segment-intron 10000 -G $genome.gene.gff3 --transcriptome-index=transcriptome.index -x 1 -p $thread  -o $tgtfolder/$sam $genome ", join(" ", @R), "\n";
		
		print SHL "touch $dir/$sam.done.log\n";
		close(SHL);
		system("chmod u+x $shell");
		system("qsub $shell");
	}elsif($platform eq "zcluster"){
		print SHL "mkdir -p $tgtfolder/$sam\n";
		print SHL "time /usr/local/tophat/2.0.13/bin/tophat -i 30 -I 10000 -g 1 -F 0.05 $q --min-segment-intron 30 --max-segment-intron 10000 -G $genome.gene.gff3 --transcriptome-index=transcriptome.index -x 1 -p $thread  -o $tgtfolder/$sam $genome ", join(" ", @R), "\n";
		
		print SHL "touch $dir/$sam.done.log\n";
		close SHL;
		system("chmod u+x $shell");
		system("qsub -q rcc-30d -pe thread $t -l mem_total=20g $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
}

system("echo 'Finished 06.folder.tophat.pl!' >> job.monitor.txt");


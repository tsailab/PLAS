#!/bin/perl -w
use strict;

my $working_dir = $ARGV[0];
my $dat_dir     = $ARGV[1];

## paired
my $dat_dir_paired = $dat_dir;
opendir(ARGV,"$dat_dir_paired");
my @files=grep(/\_CA\_R.$/,readdir ARGV);
my @file_sort = sort @files;
my $num_sample = (@file_sort+0) / 2;
my @read1 = ();
my @read2 = ();
for (my $index = 0; $index < $num_sample ;$index ++ ){	   
	   my $file1 = $file_sort[2*$index];
	   my $file2 = $file_sort[2*$index+1];	   
	   push @read1, $file1;
	   push @read2, $file2;
}
	  
### wirte script	   
open  TGT, ">run_all_trinity.sh"; 
print TGT '#!/bin/bash',"\n";
print TGT '#PBS -N j_s_trinity',"\n";
print TGT '#PBS -q batch',"\n";
print TGT '#PBS -l nodes=1:ppn=24:HIGHMEM',"\n";
print TGT '#PBS -l walltime=600:00:00',"\n";
print TGT '#PBS -l mem=110gb',"\n";
print TGT 'module load trinity/2.0.6-UGA',"\n";

print TGT 'cd ',$working_dir,"\n";
my $out   = "ShanmuEMB_trinity";
 my $left  = join(',',@read1 );
my $right  = join(',',@read2 );
	     
print TGT 'time Trinity  --output Trinity_Shanmu --CPU 12 --seqType  fq  --max_memory 110G --no_version_check'," \\","\n",
           ' --left ',$left, " \\\n",
	        ' --right ',$right, " \\\n",
	        ' 1>job.out 2>job.err  ', " \\\n";  

 __END__
 








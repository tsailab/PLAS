#!/bin/perl -w
use strict;

my $working_dir =  $ARGV[0];
my $dat_dir     =  $ARGV[1];

open  TGT_ALL, ">run_all_blat.sh"; 
print TGT_ALL '#!/bin/sh',"\n";
print TGT_ALL 'cd ',$working_dir,"\n";

opendir(ARGV,"$dat_dir");
my @files=grep(/\.fastq$/,readdir ARGV);

my $num_c =0;
foreach my $file(sort @files){	
	   $num_c++;
	   my $sh_file = "r".$num_c.'_'.$file.'.sh';	 
	
	   ## sh files
	   
	   
	   open  TGT, ">$sh_file"; 
	 	  print TGT '#!/bin/bash',"\n";
	 	  print TGT '#PBS -N j_blat',"\n";
	 	  print TGT '#PBS -q batch',"\n";
	 	  print TGT '#PBS -l nodes=1:ppn=1:HIGHMEM',"\n";
	 	  print TGT '#PBS -l walltime=24:00:00',"\n";
	 	  print TGT '#PBS -l mem=50gb',"\n";
	   


	      print TGT 'cd ',$working_dir,"\n";
	  
	      my $file_in= $dat_dir.'/'.$file;
	     
	      print TGT 'perl step_1_Fastq_vs_rRNA_only_Sapelo.pl ',$file_in,"\n";
	  
	     ### ALL input files
	       print TGT_ALL 'qsub  ', $sh_file,"\n";

 }






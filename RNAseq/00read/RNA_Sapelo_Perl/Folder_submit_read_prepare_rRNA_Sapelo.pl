#!/bin/perl -w
use strict;

my $working_dir = $ARGV[0];
my $dat_dir     = $ARGV[1];

open  TGT_ALL, ">run_all_prepare.sh"; 
print TGT_ALL '#!/bin/sh',"\n";
print TGT_ALL 'cd ',$working_dir,"\n";

opendir(ARGV,"$dat_dir");
my @files=grep(/\.fastq$/,readdir ARGV);

my @file_sort = sort @files;

my $num_sample = (@file_sort+0) / 2;


my $num_c =0;
for (my $index = 0; $index < $num_sample ;$index ++ ){
	   $num_c++; 
	   
	   my $file1 = $file_sort[2*$index];
	   my $file2 = $file_sort[2*$index+1];
	   
	   my $sh_file = "r".$num_c.'_'.$file1.'.sh';	 
	   
	   if(($file1 =~ /\_R1\_/ and $file2 =~ /\_R2\_/)  or ($file1 =~ /\_1\_/ and $file2 =~ /\_2\_/) or ($file1 =~ /\-1\_CA/ and $file2 =~ /\-2\_CA/) ){
	   	     # correct order
	   }
	   else{
	   	    die "not correct order  $file1 $file1\n";
	   }
	
	   ## sh files
	      open  TGT, ">$sh_file"; 
	 	  print TGT '#!/bin/bash',"\n";
	 	  print TGT '#PBS -N j_Read_prepare',"\n";
	 	  print TGT '#PBS -q batch',"\n";
	 	  print TGT '#PBS -l nodes=1:ppn=1:HIGHMEM',"\n";
	 	  print TGT '#PBS -l walltime=24:00:00',"\n";
	 	  print TGT '#PBS -l mem=40gb',"\n";
	 	  
	 	  
	      print TGT 'cd ',$working_dir,"\n";
	  
	      #my $file_in= $dat_dir.$file;
	     
	      print TGT 'perl step_2_Reads_Prepare_v2.pl  ',$file1,' ',$file2,' ',$dat_dir,' ', $working_dir,"\n";
	  
	   ### ALL input files
	      print TGT_ALL 'qsub ', $sh_file,"\n"; 
 }







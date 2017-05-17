#!/bin/perl -w

use strict;
### First Fastq file
my $filein1=$ARGV[0];
### Second Fastq file
my $filein2=$ARGV[1];
####
my $dirIn_1= $ARGV[2];
my $out_file = $ARGV[3];

open (TGT, ">$out_file")||die "Cannot write $out_file \n";

#######################################################################################
#------------------check the fastq version
open(SRC1,"$dirIn_1/$filein1")||die "can not open file:";

my $is_version_1_8 = '';

for(1..20){
	my $line=<SRC1>;
	if($line=~/^\@/){
		if($line =~ m/\s\d/){
			  $is_version_1_8  = 1;
		}
		elsif($line =~ m/\/\d/){
			  $is_version_1_8  = 0;
		}
	    elsif($line =~ m/\slength/){
			  $is_version_1_8  = 2;
		}
		else{
			print "New format again\n";
			die;
		}
		last;
	}
}
print "Is 1.8 version :", $is_version_1_8,"\n";

close(SRC1);

open(SRC1,"$dirIn_1/$filein1")||die "can not open file:";
open(SRC2,"$dirIn_1/$filein2")||die "can not open file:";

my $read_total_num = 0;
my $read_same_num  = 0;

my $read_1_incorrent_num = 0;
my $read_2_incorrent_num = 0;


while (my $line=<SRC1>) 
{
    if($line=~/^\@/)
    {   
    	chomp $line;
    	my $id = substr $line,1;
    	my $id_out = '';
    	if($is_version_1_8 !=0){
    	     my @id_s    = split(/\s+/, $id);
    		  $id_out  = $id_s[0];
    	}
    	else{
    		  $id_out = substr $id,0,-2;
    	}
    	
       ## We are interested in ID only ignore other lines
        $line=<SRC1>;
        $line=<SRC1>;
        $line=<SRC1>;
      
        $line=<SRC2>;
        if(eof (SRC2)){
        	print "Read2 file has less lines\n";
        	die;
        }
        chomp $line;
        
        ####
        my $id_2 = substr $line,1;
    	my $id_out_2 = '';
    	if($is_version_1_8 !=0){
    	     my @id_s    = split(/\s+/, $id_2);
    		  $id_out_2  = $id_s[0];
    	}
    	else{
    		  $id_out_2 = substr $id_2,0,-2;
    	}


        ## We are interested in ID only ignore other lines
        $line=<SRC2>;     
        if(eof (SRC2)){
        	print "Read2 file has less lines\n";
        	die;
        }
        $line=<SRC2>;     
        if(eof (SRC2)){
        	print "Read2 file has less lines\n";
        	die;
        }
        
        $line=<SRC2>;     
      

        ##### compare read1  and read2 
        
        $read_total_num ++;
        
        if($id_out eq $id_out_2){
        	   $read_same_num++;
        }
        
        if($is_version_1_8 == 1){
    	      if($id =~ m/\s2/){
    	      	   $read_1_incorrent_num++;
    	      }
    	      
    	      if($id_2 =~ m/\s1/){
    	      	   $read_2_incorrent_num++;
    	      }
    	}
    	else{
    		  if($id =~ m/\/2/){
    	      	   $read_1_incorrent_num++;
    	      }
    	      
    	      if($id_2 =~ m/\/1/){
    	      	   $read_2_incorrent_num++;
    	      }
    	}
    	
    }## @

}

print TGT "Total number \t$read_total_num\n";
print TGT "Same order read number \t$read_same_num\n";

print TGT "Incorrect reads in reads 1 \t$read_1_incorrent_num\n";
print TGT "Incorrect reads in reads 2 \t$read_2_incorrent_num\n";


if(not eof (SRC2)){
      print "Read2 file has more lines\n";  	
}

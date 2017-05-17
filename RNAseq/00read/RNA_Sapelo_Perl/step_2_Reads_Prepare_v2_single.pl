#!/bin/perl -w

use strict;
### First Fastq file
my $filein1 = $ARGV[0];
####
my $dirIn_1 = $ARGV[1];
my $dirIn_2 = $ARGV[2];

#3, directory of read1 and read2
#4, directory of the output file from Fastq_vs_Mit_Chl_rRNA.pl

###########################################################################################

#--- check file names
my $blat_hit1 = '';

my $file_out1=(substr $filein1,0,-6).'_R1';

if($filein1 =~ m/\.fastq/i){
	 $blat_hit1=(substr $filein1,0,-6).'_mrc_hit.tab';
}
elsif($filein1 =~ m/\.fq/i){
	 $blat_hit1=(substr $filein1,0,-3).'_mrc_hit.tab';
     $file_out1=(substr $filein1,0,-3).'_R1';  
}
else{
     $blat_hit1=(substr $filein1,0).'_mrc_hit.tab';
     $file_out1=(substr $filein1,0).'_R1';
}

#######################################################################################
#----------------------------------------------------
## Check blat result

my %blat_out_1=();
print 'Built blat index',"\n";

open(SRC0,"$dirIn_2/$blat_hit1")||die "can not open file:";
while (my $line=<SRC0>) 
{
    chomp $line;
    my @temp=split(/\t/,$line);
    $blat_out_1{$temp[0]}=1;
}


my $read1_filter_num = 0;

#--------------- reads index
my %read_id1=();
print 'Build read 1 index',"\n";
open(SRC1,"$dirIn_1/$filein1")||die "can not open file:";

#------------------check the fastq version
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

open(TGT1,">$file_out1")||die "can not open file:";
while (my $line=<SRC1>) 
{
    chomp $line;
    if($line=~/^\@/)
    {   
    	my $id = substr $line,1;
    	my $id_out = '';
    	if($is_version_1_8 != 0){
    	   my @id_s    = split(/\s+/, $id);
    		  $id_out  = $id_s[0];
    	}
    	else{
    		  $id_out = substr $id,0,-2;
    	}

    	
    	my @read_item=();
    	push @read_item,$line;    ###   1
    	
    	$line=<SRC1>;
    	chomp $line;
    	push @read_item,$line;    ####  2
    	
    	my $len=length($line);
    	
    	$line=<SRC1>;
    	chomp $line;
    	push @read_item,$line;    ####  3
    	
    	$line=<SRC1>;
    	chomp $line;
    	push @read_item,$line;    ####  4
    	
    	#checking blat
    	if(exists $blat_out_1{$id_out})
    	{
    		$read1_filter_num++;
    		next;
    	}
    	
    	if($len>20)
    	{
    		#$read_id1{$id_out}=\@read_item;
    		
    		print TGT1 join("\n",@read_item),"\n";
    	}
    }
}





print "Read1 filtered number\t", $read1_filter_num,"\n";


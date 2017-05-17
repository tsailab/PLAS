#!/bin/perl -w

use strict;
### First Fastq file
my $filein1=$ARGV[0];
### Second Fastq file
my $filein2=$ARGV[1];
####
my $dirIn_1 = $ARGV[2];
my $dirIn_2 = $ARGV[3];

#3, directory of read1 and read2
#4, directory of the output file from Fastq_vs_Mit_Chl_rRNA.pl

###########################################################################################

#--- check file names
my $blat_hit1 = '';
my $blat_hit2 = '';
my $file_out1=(substr $filein1,0,-6).'_R1';
my $file_out2=(substr $filein2,0,-6).'_R2';

if($filein1 =~ m/\.fastq/i){
	 $blat_hit1=(substr $filein1,0,-6).'_mrc_hit.tab';
     $blat_hit2=(substr $filein2,0,-6).'_mrc_hit.tab';
}
elsif($filein1 =~ m/\.fq/i){
	 $blat_hit1=(substr $filein1,0,-3).'_mrc_hit.tab';
     $blat_hit2=(substr $filein2,0,-3).'_mrc_hit.tab';
     $file_out1=(substr $filein1,0,-3).'_R1';
     $file_out2=(substr $filein2,0,-3).'_R2';
     
     
}
else{
     $blat_hit1=(substr $filein1,0).'_mrc_hit.tab';
     $blat_hit2=(substr $filein2,0).'_mrc_hit.tab';
     $file_out1=(substr $filein1,0).'_R1';
     $file_out2=(substr $filein2,0).'_R2';
     
}

#######################################################################################
#----------------------------------------------------
## Check blat result

my %blat_out_1=();
my %blat_out_2=();
print 'Built blat index',"\n";

open(SRC0,"$dirIn_2/$blat_hit1")||die "can not open file:";
while (my $line=<SRC0>) 
{
    chomp $line;
    my @temp=split(/\t/,$line);
    $blat_out_1{$temp[0]}=1;
}

open(SRC0,"$dirIn_2/$blat_hit2")||die "can not open file:";
while (my $line=<SRC0>) 
{
    chomp $line;
    my @temp=split(/\t/,$line);
    $blat_out_2{$temp[0]}=1;
}

my $read1_filter_num = 0;
my $read2_filter_num = 0;

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


while (my $line=<SRC1>) 
{
    chomp $line;
    if($line=~/^\@/)
    {   
    	my $id = substr $line,1;
    	my $id_out = '';
    	if($is_version_1_8 == 1){
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
    		$read_id1{$id_out}=\@read_item;
    	}
    }
}



open(TGT1,">$file_out1")||die "can not open file:";
open(TGT2,">$file_out2")||die "can not open file:";

my @low_reads=qw/
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
/;

if($is_version_1_8 == 1){
	 @low_reads= qw|
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
     +
     |;
     
     push @low_reads,'#############################################';
}

print 'Build read 2 index',"\n";
open(SRC1,"$dirIn_1/$filein2")||die "can not open file:";
while (my $line=<SRC1>) 
{
    chomp $line;
    if($line=~/^\@/)
    {   
    	my $id= substr $line,1;
    	
    	my $id_out = '';
    	if($is_version_1_8 == 1){
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
    	if(exists $blat_out_2{$id_out})
    	{
    		$read2_filter_num++;
    		next;
    	}
    	
    	if($len>20)
    	{
  
    		#---  output---
    		
    		if(exists $read_id1{$id_out})
    		{    ## paired
    			my @read_1=@{$read_id1{$id_out}};
    			print TGT1 join("\n",@read_1),"\n";
    			print TGT2 join("\n",@read_item),"\n";
    			delete($read_id1{$id_out});
    		}
    		else{
    			## read2 only 
    			## get the fake read 1 name
    			 my $read1_name = $id;
    		     if($is_version_1_8 == 1){
    	              if(not $read1_name =~ s/\s2/ 1/g){
    	              	   print "1.8 read2 to read1 error\n";
    	              	   die;
    	              }
    	         }
    	         else{
    		          if(not $read1_name =~ s/\/2/\/1/g){
    		          	   print "1.5 read2 to read1 error\n";
    		          	   die;
    		          }
    	         }
    			
    			print TGT1 '@',$read1_name,"\n";
    			print TGT1 join("\n",@low_reads),"\n";
    			print TGT2 join("\n",@read_item),"\n";
    			
    		}
    	}
    }
}

### read 1 only

for my $id(keys %read_id1){
	  my @read_1=@{$read_id1{$id}};
	  print TGT1 join("\n",@read_1),"\n";
	 ## id1 to id2
      
      my $read2_name = $read_1[0];
      
      if($is_version_1_8 == 1){
    	      if(not $read2_name =~ s/\s1/ 2/g){
    	      	   print "1.8 read1 to read2 error\n";
    	      	   die;
    	      }
      }
      else{
    		  if(not $read2_name =~ s/\/1/\/2/g){
    		  	   print "1.8 read1 to read2 error\n";
    		  	   die;
    		  }
      }
	 

    print TGT2 $read2_name,"\n";
    print TGT2 join("\n",@low_reads),"\n";
	
}

print "Read1 filtered number\t", $read1_filter_num,"\n";
print "Read2 filtered number\t", $read2_filter_num,"\n";


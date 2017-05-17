#!/bin/perl -w

use strict;
####  First input file: FastQ file
my $file_full   =    $ARGV[0];
my @ooo = split(/\//,$file_full);
my $filein   =    $ooo[-1];

###  settings 
my $blat_db = '/home/lxue/database/rRNA_only_NR.fas';
#my $ooc_file = '/panfs/pstor.storage/escratch1/lxue_Jun_01/SUT_DF/02mit_rrna/11.ooc';
my $commad = '/usr/local/apps/blat/latest/bin/blat';


print 'The blat software is here: ',$commad,"\n";
print 'The blat database is here: ',$blat_db,"\n";

###########################################################
#--- check file names
my $file_fas = '';
if($filein =~ m/\.fastq/i){
	$file_fas  = (substr $filein,0,-5).'fasta';
}
elsif($filein =~ m/\.fq/i){
	$file_fas  = (substr $filein,0,-2).'fasta';
}
else{
	$file_fas  = $filein.'.fasta';
}

####################################################################
#-------------  fastq  to fasta  --------------------

open(SRC,"$file_full") || die;
open(TGT,">$file_fas") || die "can not WRITE file:";


###### check fastq version
my $is_version_1_8 = '';
for(1..20){
	my $line=<SRC>;
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
close(SRC);

################################
## Re open the file

open(SRC,"$file_full") || die;
	
while (my $line=<SRC>)
{
    chomp $line;
    if($line=~/^\@/)
    {
    	my $id     = substr $line,1;
    	my $id_out = '';
    	
    	if($is_version_1_8  == 1){
    	   my @id_s    = split(/\s+/, $id);
    		  $id_out  = $id_s[0];
    	}
    	elsif($is_version_1_8 == 2){
    	   my @id_s    = split(/\s+/, $id);
    		  $id_out  = $id_s[0];
    	}
    	else{
    		  $id_out = substr $id,0,-2;
    	}

    	$line=<SRC>;
    	chomp $line;

    	my $len=length($line);
    
    	if($len>20)
    	{
    	   print TGT ">", $id_out,"\n";
    	   print TGT "$line\n";
    	}
    	$line=<SRC>;    		
    	$line=<SRC>;
    }
 }
 
 close(TGT);
 
 ##############################################################################################
 #--------- run  blat ----------------------------
 my $blat_in = $file_fas;

 my $blat_out = (substr $blat_in,0,-6).'_mrc.psl';

 `$commad  $blat_db $blat_in  $blat_out `;
 
 #  Ptri_156_cds_len.fas Ptri_156_cds_len_p4.fas -ooc= Ptri_156_cds_len_p4_sel.bls -out=blast8
 
 
################################################################################################    
 #---------------check blat result
   
my $sourcefile = $blat_out;
my $targetfile = substr ($blat_out,0,-4).'_hit.tab';
   
open(SRC,"$sourcefile") || die "can not open file:".$sourcefile;
open(TGT,">$targetfile") || die "can not write file:".$targetfile;

my %blast_first=();
my $cut_off = 0.9;

## remove first 5 lines
for(1..5){
	my $line=<SRC>;
	#print OUT $line;
}

while (my $line=<SRC>)
{
     chomp $line;
     my @lineArr = split(/\t/,$line);
     
     	#my $scaffold = $lineArr[13];
		my $read_id  = $lineArr[9];
		my $q_size   = $lineArr[10];
		my $q_match  = $lineArr[0];
		
		my $target  = $lineArr[13];
		
		my $ratio  = $q_match/ $q_size;
		if( $ratio < $cut_off){
			next;
		}
		
     	if(not exists $blast_first{$read_id})
     	{
     		 $blast_first{$read_id} = 1;
     		 print TGT $read_id,"\t",$ratio,"\t",$target,"\n";
     	}
}

close(SRC);
close(TGT);

`rm $blat_out`;
`rm $blat_in`;

   
   
   
   
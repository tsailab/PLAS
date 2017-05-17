#!/bin/perl -w

#input1 input file 
#input2 the output file
use strict;
my $src    =$ARGV[0];
my $fileout=$ARGV[1];

open(SRC,"$src")||die "can not open file:";
open(TGT,">$fileout")||die "can not open file:";

my $hash = ();
while(<SRC>) 
{
	if(m/^#/)
	{
		next;
	}
	
    my @line = split(/\t/,$_);
    
    if($line[8] =~ m/([^\r\n]+)/){
    	  $line[8] = $1;
    }    
    my $id='';
    my $trans = '';
    my $gene = '';
       if($line[8] =~m/ID\=([^\;]+)/){
    		$id=$1;
    	}
    	if($line[8] =~m/Genbank\:([^\;]+)/){
    		$trans =$1;
    	}
    	if($line[8] =~m/gene\=([^\;]+)/){
    		$gene =$1;
    	}
    
    if(uc($line[2]) eq 'CDS'){	
    	  pop @line;
		  
		  if($gene ne ""){
		  	if (not exists $$hash{$gene}){
		  		 print TGT $gene,"\t",$trans,"\n";
		  		  $$hash{$gene}= 1;
		  	}
		    
		  }
    }
} ## while file

print "Read GFF done\n";
    	  
#!/bin/perl -w

use strict;

my $srcfolder = $ARGV[0];
my $targetfile= $ARGV[1];
open(TGT,">$targetfile")||die "can not WRITE file:";
opendir(ARGV,"$srcfolder");
my @files=grep(/\.tab$/,readdir ARGV);

my %hit_hash = ();
my %hit_names = ();
foreach my $file(@files){	
	open (SRC,"$srcfolder\/$file");
   #each file
    # <SRC>;
      my $chr=substr $file,0,-4;

      while (my $name=<SRC>) { 
      	  chomp $name;  
      	  my @temp = split(/\t/,$name);
      	  my $hit  = uc($temp[2]);
      	  if($hit =~ m/MARKER/){
      	  	   $hit_hash{$chr}{$hit}{$temp[0]} = 1;
      	  	   $hit_names{$hit} = 1
      	  }
      	  else{
      	  	  $hit_hash{$chr}{"rRNA"}{$temp[0]} = 1;
      	  	  $hit_names{"rRNA"} = 1
      	  }
     }
     

}

my @hit_names = sort keys %hit_names;
print TGT "file\t",join("\t",@hit_names),"\n";

for my $file (keys %hit_hash ){
	 my @out = ();
	for my $hit (@hit_names){
		my $num = 0;
		if(exists  $hit_hash{$file}{$hit}){
			my @ooo = keys %{$hit_hash{$file}{$hit}};
			$num = @ooo +0;
		}
		push @out,$num;
	}
	print TGT $file,"\t",join("\t",@out),"\n";
}
    

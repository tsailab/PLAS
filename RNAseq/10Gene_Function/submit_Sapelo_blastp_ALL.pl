#!/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Seq;
my $out;
#in a folder
my $working =  $ARGV[0];
my $srcfile  =  $ARGV[1];
my $each = 2000;
my $split_out = $ARGV[2];
my $database =  $ARGV[3];

     
my $seqio_object = Bio::SeqIO->new(-file   => "$srcfile", -format => "fasta");  

my $i=0;
my $file;
my $oldfile=0;
my $seqout_object; 
while (my $seq_object = $seqio_object->next_seq){
    $i++;   
    $file=int(($i-1)/$each)+1;
    if($file != $oldfile){
    	my $outfile=$split_out.$file.'.fas';
    	open TGT, ">$outfile"; 
     #$seqout_object = Bio::SeqIO->new(-file   => ">$outfile", -format => "fasta"); 
    	$oldfile=$file;
    	print TGT ">",$seq_object->id,"\n",$seq_object->seq,"\n";
     #$seqout_object->write_seq($seq_object);
    }    
   else{
      print TGT ">",$seq_object->id,"\n",$seq_object->seq,"\n";
     #$seqout_object->write_seq($seq_object);
   }  
}

my $max=$file;

    

open  TGT_ALL, ">run_blast".$split_out.".sh"; 
print TGT_ALL '#!/bin/sh',"\n";
print TGT_ALL 'cd ',$working,"\n";


for(my $idx=1;$idx<=$max;$idx++){ 
	 my $fas_file=$split_out.$idx.'.fas';
	
	

	 	  my $sh_file="r".$idx.'_'.$split_out .'.sh';
	 	  my $bls_out = $split_out.$idx.'_prot'.'.bls';
	 	  print TGT_ALL 'qsub ', $sh_file,"\n"; 
	 	  
	       # shell file
	 	  open  TGT, ">$sh_file"; 
	 	  print TGT '#!/bin/bash',"\n";
	 	  print TGT '#PBS -N j_blastp'.$idx,"\n";
	 	  print TGT '#PBS -q batch',"\n";
	 	  print TGT '#PBS -l nodes=1:ppn=5:HIGHMEM',"\n";
	 	  print TGT '#PBS -l walltime=24:00:00',"\n";
	 	  print TGT '#PBS -l mem=20gb',"\n";
	 	  
	 	  print TGT 'module load  ncbiblast+/2.2.29',"\n";
	 	  
	 	  print TGT 'cd ',$working,"\n";
	 	  
	 	  print TGT 'time blastp -query ', 
	 	  $fas_file, ' -db ',$database,' -out ', $bls_out,' -evalue 1e-05 -num_threads 5   -num_descriptions 1 -num_alignments 1',"\n";
	 	
}


	














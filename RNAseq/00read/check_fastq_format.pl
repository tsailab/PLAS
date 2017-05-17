#!/bin/perl -w
use strict;
use List::MoreUtils qw( minmax );
my $dat_dir  = $ARGV[0];
my $out_file     = $ARGV[1];

open TGT, ">$out_file" || "die cann't write file\n";
opendir(ARGV,"$dat_dir");
my @files=grep(/\.fastq$|\.fq$/,readdir ARGV);

my @file_sort = sort @files;


my $num_c =0;
for my $file1 (@file_sort) {
	   $num_c++; 
	   ####
	   my $read1 = $dat_dir .'/'.$file1;
       my ($qual_format,$length_format) = &get_read_format($read1);
       
       print TGT  $file1,"\t",$qual_format,"\t",$length_format,"\n";
}

          
	                
	                
sub	  get_read_format{
	my $read_in  = shift @_;
	my $limit = 200;
    open FILEIN, "<$read_in";
    my ($min,$max);
    my $cnt =0;
    my $len_max = 0;
	while (my $id = <FILEIN>) {
			$id =~ m/^@/ || die "expected @ not found in line 1!\n";
			my $seq = <FILEIN>;
			my $sep = <FILEIN>;
			$sep =~ m/^\+/ || die "expected + not found in line 3!\n";
			my $qual = <FILEIN>;
			chomp($qual);
			if(length($qual)>$len_max){
				$len_max = length($qual);
			}
			$cnt++;
			$cnt>=$limit && last;
			# char to ascii
			my @chars = split("", $qual);
			my @nums = sort { $a <=> $b } (map { unpack("C*", $_ )} @chars);
			if ($cnt==1) {
			     ($min, $max) = minmax @nums;
			} else {
					my ($lmin, $lmax) = minmax @nums; # local values for this read
					$lmin<$min ? $min=$lmin : $min=$min;
					$lmax>$max ? $max=$lmax : $max=$max;
			}
     }# while
     
     my $qual = 33;
     if ($min<33 || $max>105) { die "Quality values corrupt. found [$min; $max] where [33; 104] was expected\n"; }

	 if ($min >=59 && $max <= 110 )  {
			$qual = 64;
	 }
	 if ($min>=33 && $max <= 74 ) {
	 	$qual = 33;
	 }
	 return ($qual,$len_max);
}



__END__
    	   # 1.8 version fastq   
	      print TGT  '/usr/local/tophat/2.0.4/bin/tophat  -i 30 -I 10000 -g 50 -F 0 --min-segment-intron 30 --max-segment-intron 10000  \\',"\n",
	          '-G  /panfs/pstor.storage/grphomes/cjtlab/lxue/genome/Ptrichocarpa_v2.2mMay.gff3   -p 3 \\',"\n",
	          '-o ',$out,' --bowtie1  \\',"\n",
	          ' /panfs/pstor.storage/grphomes/cjtlab/lxue/genome/Ptrichocarpa_156_bowtie \\',"\n",
	          $read1,' ',$read2,"\n"; 
	  
               

##!/bin/bash
#cd /panfs/pstor.storage/escratch/lxue_Aug_03/10cufflinks/01reads
#/usr/local/tophat/latest/bin/tophat  -i 30 -I 10000 -g 50 -F 0.05  --min-segment-intron 30 --max-segment-intron 10000  \
#-G  /panfs/pstor.storage/escratch/lxue_Aug_03/10v3_genome/Ptrichocarpa_v3.0_210_gene.gff3  -p 4 \
#-o  T3_SH1 \
#/panfs/pstor.storage/escratch/lxue_Aug_03/10v3_genome/Ptrichocarpa_v3.0_210 \
#T3_SH1a_1_1_CA4_R1,T3_SH1a_1_2_CA4_R1,T3_SH1a_1_3_CA4_R1,T3_SH1a_1_4_CA4_R1 \
#T3_SH1a_2_1_CA4_R2,T3_SH1a_2_2_CA4_R2,T3_SH1a_2_3_CA4_R2,T3_SH1a_2_4_CA4_R2 ;




cd /panfs/pstor.storage/escratch/lxue_Aug_03/10cufflinks/01reads
/usr/local/tophat/latest/bin/tophat  -i 30 -I 10000 -g 50 -F 0.05  --min-segment-intron 30 --max-segment-intron 10000  \
-G  /panfs/pstor.storage/escratch/lxue_Aug_03/05v2.2_genome/Ptrichocarpa_156_gene.gff3  -p 4 \
-o  T3_SH1 \
/panfs/pstor.storage/escratch/lxue_Aug_03/05v2.2_genome/Ptrichocarpa_156 \
T3_SH1a_1_1_CA4_R1,T3_SH1a_1_2_CA4_R1,T3_SH1a_1_3_CA4_R1,T3_SH1a_1_4_CA4_R1 \
T3_SH1a_2_1_CA4_R2,T3_SH1a_2_2_CA4_R2,T3_SH1a_2_3_CA4_R2,T3_SH1a_2_4_CA4_R2 ;








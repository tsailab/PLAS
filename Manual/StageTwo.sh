#!/bin/bash

platform='Zcluster'
mode='paired-end'

#time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 4 -max_target_seqs 1
#wait
#time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.v1.fasta -dbtype nucl
#wait
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.v1.fasta -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5
wait
time perl 00.script/c11.remove.redundancy.pl 08.full.length/Final.v1.blastn.xml.out 08.full.length/Final.v1.fasta 08.full.length/Final.v2.fasta self 08.full.length/Final.v1.ptr.blastx.out > 08.full.length/record.run1
wait
time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.v2.fasta -dbtype nucl
wait
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.v2.fasta -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5
wait
time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 4 -max_target_seqs 1
wait
time perl 00.script/101.transfer.saturate.seq.pl 08.full.length/Final.v2.ptr.blastx.out 08.full.length/Final.v2.fasta 01.data/00.PriorData/ptr.proteome.fa 08.full.length Final.v2.ptr pct 0.02
wait
###MADE IT TO HERE
cd 08.full.length/
ln -sf Final.v2.ptr.full.contigs.nucl.fasta Final.fasta
cd ../

#WORK OUT PERL JOB SUBMISSION TIMING AS WELL#
########
#module load bowtie2/2.2.4
#export PATH=$PATH:/usr/local/bowtie2/2.2.3/bin/
#time bowtie2-build -f -q 08.full.length/Final.fasta 08.full.length/Final
#wait
#time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 08.full.length/Final 09.bowtie.full.length unmap $platform
#wait
#time perl 00.script/c10.unmapped.reads.trinity.pl 09.bowtie.full.length $platform
#wait
#time perl 00.script/06.truncate.header.pl 10.unmapped.reads.trinity/Trinity.fasta 10.unmapped.reads.trinity/Trinity.new.fasta
########
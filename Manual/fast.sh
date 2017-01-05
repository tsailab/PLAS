#!/bin/bash

platform='Zcluster'
mode='paired-end'

export PATH=$PATH:/usr/local/bowtie2/2.2.3/bin/
time bowtie2-build -f -q 08.full.length/Final.fasta 08.full.length/Final
wait
time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 08.full.length/Final 09.bowtie.full.length unmap $platform
wait
time perl 00.script/c10.unmapped.reads.trinity.pl 09.bowtie.full.length $platform

while [ ! -f flag.txt ]
do
   sleep 2
done

time perl 00.script/06.truncate.header.pl 10.unmapped.reads.trinity/Trinity.fasta 10.unmapped.reads.trinity/Trinity.new.fasta

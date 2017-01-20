#!/bin/bash
export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
export PATH=/usr/local/gmap-gsnap/latest/bin/:${PATH}

time /usr/local/trinity/r20140717/Trinity --seqType fa --CPU 4 --JM 20G --left 09.bowtie.full.length/1T/unmapped.reads.1T.1.fasta,09.bowtie.full.length/2T/unmapped.reads.2T.1.fasta,09.bowtie.full.length/3T/unmapped.reads.3T.1.fasta --right 09.bowtie.full.length/1T/unmapped.reads.1T.2.fasta,09.bowtie.full.length/2T/unmapped.reads.2T.2.fasta,09.bowtie.full.length/3T/unmapped.reads.3T.2.fasta --output 10.unmapped.reads.trinity
echo You may proceed! >> FLAGFILE.txt
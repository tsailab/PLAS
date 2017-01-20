#!/bin/bash
export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
export PATH=/usr/local/gmap-gsnap/latest/bin/:${PATH}

time /usr/local/trinity/r20140717/Trinity --seqType fa --CPU 8 --JM 400G --left 01.data/02.Fasta/1T/1T.combined.R1.fasta,01.data/02.Fasta/2T/2T.combined.R1.fasta,01.data/02.Fasta/3T/3T.combined.R1.fasta --right 01.data/02.Fasta/1T/1T.combined.R2.fasta,01.data/02.Fasta/2T/2T.combined.R2.fasta,01.data/02.Fasta/3T/3T.combined.R2.fasta --output 05.Trinity

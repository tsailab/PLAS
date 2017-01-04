#!/bin/bash

platform='Zcluster'
mode='paired-end'
echo "Process 1 begun"
time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 4 -max_target_seqs 1
wait
echo "Process 1 complete!"
echo "Process 2 begun!"
time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.v1.fasta -dbtype nucl
echo "Process 2 complete!"

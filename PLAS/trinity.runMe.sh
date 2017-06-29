#!/bin/bash
run="$1"
sub="$2"
mode="paired-end"
thread="2"
memory="24"
tgtfolder=06.assembly/03.bowtie.nucl/run."$run"/"$sub".trinity
srcfolder=04.retrieve.reads/03.bowtie.nucl/run."$run"/"$sub"
mkdir $tgtfolder
module load trinityrnaseq/2.2.0
module load java

time Trinity --seqType fa --CPU $thread --max_memory "$memory"G \
        --left "$srcfolder"/retrieved."$sub".R1.fasta --right "$srcfolder"/retrieved."$sub".R2.fasta \
        --output $tgtfolder --min_contig_length 100
		

echo "Trinity output to $tgtfolder/$sub.trinity"

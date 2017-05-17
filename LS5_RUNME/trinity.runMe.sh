#!/bin/bash
b="$1"
sub="$2"
srcfolder="04.retrieve.reads/03.bowtie.nucl/run.$b"
tgtfolder="06.assembly/03.bowtie.nucl/run.$b"
reffolder="01.data/05.SplitGenes/01.Protein/run.0"
scale="genome"
mode="paired-end"
thread="1"
memory="1"
module load bowtie/1.1.2
echo "Memory before Trinity: " >> job.monitor.txt
grep MemTotal /proc/meminfo | awk '{print} $2' >> "job.monitor.txt"
if [ $mode == "paired-end" ]; then
	time Trinity --seqType fa --CPU $thread --max_memory "$memory"G --left $srcfolder/$sub/retrieved.$sub.R1.fasta --right $srcfolder/$sub/retrieved.$sub.R2.fasta --output $tgtfolder/$sub.trinity --min_contig_length 100
	echo "Memory after Trinity: " >> job.monitor.txt	
	grep MemTotal /proc/meminfo | awk '{print} $2' >> "job.monitor.txt"
	echo "Trinity output to $tgtfolder/$sub.trinity"
elif [ $mode == "single-end" ]; then
	time Trinity --seqType fa --CPU $thread --max_memory "$memory"G --single $srcfolder/$sub/retrieved.$sub.fasta --output $tgtfolder/$sub.trinity --min_contig_length 100
else
	echo "Error: No read mode available!" >> job.monitor.txt
	exit
fi

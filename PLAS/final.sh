#PBS -S /bin/bash
#PBS -q batch
#PBS -N Master_PLAS
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l walltime=12:00:00
#PBS -l mem=50gb
cd $PBS_O_WORKDIR

module load ncbiblast+
module load perl/5.20.2-thread
module load trinity/r20140717
samplesLeft=""
samplesRight=""

## Combine full length contigs from local assembly and unmapped reads assembly
for sub in 01.data/02.fasta/*; do
	if [ -d $sub ]; then
		cat 09.bowtie.full.length/${sub}/unmapped.reads.${sub}.single.fasta >> 09.bowtie.full.length/${sub}/unmapped.reads.${sub}.1.fasta
	fi
done

## Assemble unmapped reads using Trinity
## Get filenames from each sample folder

for sub in 09.bowtie.full.length/*; do 
	if [ -d $sub ]; then 
	sub=$(basename $sub)
samplesLeft+="09.bowtie.full.length/${sub}/unmapped.reads.${sub}.1.fasta,"
samplesRight+="09.bowtie.full.length/${sub}/unmapped.reads.${sub}.2.fasta,"
	fi
done

samplesLeft=${samplesLeft%?}
samplesRight=${samplesRight%?}

Trinity --seqType fa --CPU 4 --JM "20"G --left "$samplesLeft" --right "$samplesRight" --output 10.unmapped.reads.trinity

perl 00.script/06.truncate.header.pl 10.unmapped.reads.trinity/Trinity.fasta 10.unmapped.reads.trinity/Trinity.new.fasta

makeblastdb -in 08.full.length/Final.fasta -dbtype nucl

blastn -db 08.full.length/Final.fasta -query 10.unmapped.reads.trinity/Trinity.new.fasta -out 10.unmapped.reads.trinity/Trinity.new.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 1

perl 00.script/c11.remove.redundancy.pl 10.unmapped.reads.trinity/Trinity.new.blastn.xml.out 10.unmapped.reads.trinity/Trinity.new.fasta 10.unmapped.reads.trinity/temp.fasta query > 10.unmapped.reads.trinity/remove.redundancy.log

cat 08.full.length/Final.fasta 10.unmapped.reads.trinity/temp.fasta > 01.data/06.TargetTranscriptome/transcriptome.v1.fa

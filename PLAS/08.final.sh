#PBS -S /bin/bash
#PBS -N PLAS_final
#PBS -q highmem_q
#PBS -l nodes=1:ppn=12
#PBS -l mem=48gb
#PBS -l walltime=18:00:00
#PBS -d ./

cd $PBS_O_WORKDIR

#Load required modules
module load BLAST+/2.6.0-foss-2016b-Python-2.7.14
module load Perl/5.20.3-foss-2016b
module load Trinity/2.6.6-foss-2016b

samplesLeft=""
samplesRight=""

## Combine full length contigs from local assembly and unmapped reads assembly
for sub in 01.data/02.fasta/*; do
	if [ -d $sub ]; then
		cat 09.fulllength.bowtie/${sub}/unmapped.reads.${sub}.single.fasta >> 09.fulllength.bowtie/${sub}/unmapped.reads.${sub}.1.fasta
	fi
done

## Assemble unmapped reads using Trinity
## Get filenames from each sample folder

for sub in 09.fulllength.bowtie/*; do 
	if [ -d $sub ]; then 
	     sub=$(basename $sub)
         samplesLeft+="09.fulllength.bowtie/${sub}/unmapped.reads.${sub}.1.fasta,"
         samplesRight+="09.fulllength.bowtie/${sub}/unmapped.reads.${sub}.2.fasta,"
	fi
done

samplesLeft=${samplesLeft%?}
samplesRight=${samplesRight%?}

Trinity  --seqType fa --CPU 8 --max_memory "20"G --left "$samplesLeft" --right "$samplesRight" --output 10.unmappedreads.trinity

perl 00.script/06.truncate.header.pl 10.unmappedreads.trinity/Trinity.fasta 10.unmappedreads.trinity/Trinity.new.fasta

makeblastdb -in 08.full.length/Final.fasta -dbtype nucl

blastn -db 08.full.length/Final.fasta -query 10.unmappedreads.trinity/Trinity.new.fasta -out 10.unmappedreads.trinity/Trinity.new.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 1

perl 00.script/c11.remove.redundancy.pl 10.unmappedreads.trinity/Trinity.new.blastn.xml.out 10.unmappedreads.trinity/Trinity.new.fasta 10.unmappedreads.trinity/temp.fasta query > 10.unmappedreads.trinity/remove.redundancy.log

cat 08.full.length/Final.fasta 10.unmappedreads.trinity/temp.fasta > 01.data/06.TargetTranscriptome/transcriptome.v1.fa

## Script ends


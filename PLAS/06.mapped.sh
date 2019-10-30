#PBS -S /bin/bash
#PBS -N PLAS_mapped
#PBS -q highmem_q
#PBS -l nodes=1:ppn=12
#PBS -l mem=48gb
#PBS -l walltime=10:00:00
#PBS -d ./

cd $PBS_O_WORKDIR

#Load required modules
module load Perl/5.20.3-foss-2016b
module load Bowtie2/2.3.3-foss-2016b
module load BLAST+/2.6.0-foss-2016b-Python-2.7.14
module load BioPerl/1.7.1-foss-2016b-Perl-5.24.1

##Summarize all runs
# Retrieve header metadata from full length contig list
grep ">" 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta > 01.data/05.SplitGenes/03.Full.Length/count1

# Remove ">" from entries
sed -i "s/>//" 01.data/05.SplitGenes/03.Full.Length/count1

# Reformat contigs, remove redundant entries
perl 00.script/b3.full.length.format.pl 01.data/04.GeneOfInterest/GeneID.v1.txt 01.data/05.SplitGenes/03.Full.Length/count1 01.data/05.SplitGenes/03.Full.Length/count2 01.data/05.SplitGenes/03.Full.Length/count3

# assemble unmapped reads
cp 01.data/05.SplitGenes/03.Full.Length/count3 08.full.length/
cp 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta 08.full.length/
time perl 00.script/c9.get.full.length.seq.pl 08.full.length/count3 08.full.length/full.length.contigs.nucl.fasta 08.full.length/Final.v1.fasta

makeblastdb -in 01.data/00.PriorData/ptr.proteome.fa -dbtype prot
blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 32 -max_target_seqs 1

makeblastdb -in 08.full.length/Final.v1.fasta -dbtype nucl
blastn -db 08.full.length/Final.v1.fasta -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5

time perl 00.script/c11.remove.redundancy.pl 08.full.length/Final.v1.blastn.xml.out 08.full.length/Final.v1.fasta 08.full.length/Final.v2.fasta self 08.full.length/Final.v1.ptr.blastx.out > 08.full.length/record.run1

makeblastdb -in 08.full.length/Final.v2.fasta -dbtype nucl

blastn -db 08.full.length/Final.v2.fasta -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5
blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 32 -max_target_seqs 1

time perl 00.script/101.transfer.saturate.seq.pl 08.full.length/Final.v2.ptr.blastx.out 08.full.length/Final.v2.fasta 01.data/00.PriorData/ptr.proteome.fa 08.full.length Final.v2.ptr pct 0.02

cd 08.full.length/
ln -sf Final.v2.ptr.full.contigs.nucl.fasta Final.fasta
cd ../

time bowtie2-build -f -q 08.full.length/Final.fasta 08.full.length/Final

## Script ends




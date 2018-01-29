#PBS -S /bin/bash
#PBS -q batch
#PBS -N bowtie.full.length.${sub}.sh
#PBS -l nodes=1:ppn=$thread:AMD
#PBS -l walltime=48:00:00
#PBS -l mem=40gb
cd $PBS_O_WORKDIR
				
R1="01.data/02.Fasta/${sub}/${sub}.R1.fasta_pairs_R1.fasta"
R2="01.data/02.Fasta/${sub}/${sub}.R2.fasta_pairs_R2.fasta"
R3="01.data/02.Fasta/${sub}/${sub}.R1.fasta_singles.fasta"
R4="01.data/02.Fasta/${sub}/${sub}.fna"

time bowtie2 -f -x 08.full.length/Final -p 4 -1 $R1 -2 $R2 -U $R3 --un-conc 09.bowtie.full.length/${sub}/unmapped.reads.${sub}.fasta --un 09.bowtie.full.length/${sub}/unmapped.reads.${sub}.single.fasta


#PBS -S /bin/bash
#PBS -N bowtie.full.length.${sub}.sh
#PBS -q highmem_q
#PBS -l nodes=1:ppn=8
#PBS -l mem=40gb
#PBS -l walltime=08:00:00
#PBS -d ./


cd $PBS_O_WORKDIR

#Load required modules
module load Bowtie2/2.3.3-foss-2016b

sub=$SUB	
R1="01.data/02.Fasta/${sub}/${sub}.R1.fasta_pairs_R1.fasta"
R2="01.data/02.Fasta/${sub}/${sub}.R2.fasta_pairs_R2.fasta"
R3="01.data/02.Fasta/${sub}/${sub}.R1.fasta_singles.fasta"
R4="01.data/02.Fasta/${sub}/${sub}.fna"

mkdir -p 01.data/02.Fasta/$sub

time bowtie2 -f -x 08.full.length/Final -p 4 -1 $R1 -2 $R2 -U $R3 --un-conc 09.fulllength.bowtie/${sub}/unmapped.reads.${sub}.fasta --un 09.fulllength.bowtie/${sub}/unmapped.reads.${sub}.single.fasta

## Script ends

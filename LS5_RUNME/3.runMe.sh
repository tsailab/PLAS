#!/bin/bash
#SBATCH -J PLAS_runMe3
#SBATCH -o PLAS_runMe3.o%j
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p normal
#SBATCH -t 02:00:00

#SBATCH --mail-user=sjq28742@uga.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

##Additional work; check mode based on input files

#########################################################
#########derived from 01.fastq2fasta.folder.pl###########
#########################################################

echo 'Running 01.fastq2fasta.folder.pl ....' >> job.monitor.txt
mode="paired-end"
srcfolder="01.data/01.Fastq"
tgtfolder="01.data/02.Fasta"

thread="1"

rm -rf 00.script/shell.script.previous
mv 00.script/shell.script 00.script/shell.script.previous
mkdir -p 00.script/shell.script

for sub in 01.data/01.Fastq/*; do
  if [ -d "$sub" ]; then
     sub=$(basename ${sub})   
     mkdir -p 01.data/02.Fasta/$sub
        if [ $mode == "paired-end" ]; then
                awk '1 == (NR) % 4 || 2 == (NR) % 4' 01.data/01.Fastq/$sub/$sub.R1.fastq_pairs_R1.fastq | awk '{gsub("^@", ">", $0); print $0}' > 01.data/02.Fasta/$sub/$sub.R1.fasta_pairs_R1.fasta
                awk '1 == (NR) % 4 || 2 == (NR) % 4' 01.data/01.Fastq/$sub/$sub.R2.fastq_pairs_R2.fastq | awk '{gsub("^@", ">", $0); print $0}' > 01.data/02.Fasta/$sub/$sub.R2.fasta_pairs_R2.fasta
                awk '1 == (NR) % 4 || 2 == (NR) % 4' 01.data/01.Fastq/$sub/$sub.R1.fastq_singles.fastq | awk '{gsub("^@", ">", $0); print $0}' > 01.data/02.Fasta/$sub/$sub.R1.fasta_singles.fasta
                sed -i "s/\/1//g" 01.data/02.Fasta/$sub/$sub.R1.fasta_singles.fasta
                sed -i "s/\/2//g" 01.data/02.Fasta/$sub/$sub.R1.fasta_singles.fasta
        elif [ $mode == "single-end" ]; then
                awk '1 == (NR) % 4 || 2 == (NR) % 4' 01.data/01.Fastq/$sub/$sub.fastq | awk '{gsub("^@", ">", $0); print $0}' > 01.data/02.Fasta/$sub/$sub.fasta
        else
                echo 'Error: No read mode available!'
                exit 1
        fi
  fi
done

echo 'Finished converting fastq to fasta!' >> job.monitor.txt
chmod 755 -R 01.data/01.Fastq
chmod 755 -R 01.data/02.Fasta
chmod 755 -R 00.script


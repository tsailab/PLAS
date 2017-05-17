#!/bin/bash
#sBATCH -J PLAS_runMe2
#SBATCH -o PLAS_runMe2.o%j
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p normal
#SBATCH -t 01:00:00

#SBATCH --mail-user=sjq28742@uga.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

##Additional work; check mode based on input files
mode="paired-end"

#########################################################
#####derived from 01.folder.fastaCombinePairedEnd.pl#####
#########################################################

echo 'Separating paired-end and single reads...' >> job.monitor.txt
module load python

## start running the script
mkdir -p 00.script/shell.script

##Get and loop through subdirectories
for sub in 01.data/01.Fastq/*; do
  if [ -d $sub ]; then
	sub=$(basename ${sub})
        echo 'Running paired end script for $sub' >> job.monitor.txt
        time python2.7 00.script/01.fastaCombinePairedEnd.py 01.data/01.Fastq/$sub/$sub.R1.fastq 01.data/01.Fastq/$sub/$sub.R2.fastq " "
  fi
done

echo 'Finished separating reads!' >> job.monitor.txt
chmod 755 -R 01.data/01.Fastq   ##May be unnecessary
chmod 755 -R 01.data/02.Fasta
chmod 755 -R 00.script

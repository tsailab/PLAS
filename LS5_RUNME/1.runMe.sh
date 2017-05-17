#!/bin/bash
#SBATCH -J PLAS_runMe1
#SBATCH -o PLAS_runMe1.o%j
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal
#SBATCH -t 00:15:00

#SBATCH --mail-user=sjq28742@uga.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


##Additional work; check mode based on input files
mode="paired-end"
#########################################################
#######derived from 02.makeblastdb.folder.pl#############
#########################################################

echo 'Loading BLAST module...' >> job.monitor.txt
module load blast

echo 'Creating databases...' >> job.monitor.txt

##Get and loop through subdirectories

for sub in 01.data/05.SplitGenes/01.Protein/run.0/*; do
    if [ -d ${sub} ]; then # if not a directory, skip
        sub=$(basename ${sub})
        makeblastdb -in 01.data/05.SplitGenes/01.Protein/run.0/$sub/$sub.fasta -dbtype prot    
        echo "Created blast prot db in $sub..." >> job.monitor.txt

        #Make DIAMOND prot db
        $WORK/bin/diamond makedb --in 01.data/05.SplitGenes/01.Protein/run.0/$sub/$sub.fasta -d 01.data/05.SplitGenes/01.Protein/run.0/$sub/$sub
        echo "Created diamond db in $sub..." >> job.monitor.txt
    fi
done

##Get and loop through subdirectories

for sub in 01.data/05.SplitGenes/02.Transcript/run.0/*; do
    if [ -d ${sub} ]; then # if not a directory, skip
        sub=$(basename ${sub})
        makeblastdb -in 01.data/05.SplitGenes/02.Transcript/run.0/$sub/$sub.fasta -dbtype nucl
        echo "Created blast nucl db in $sub..." >> job.monitor.txt
    fi
done

echo 'Finished creating databases!' >> job.monitor.txt
chmod 755 -R 00.script
chmod 755 -R 01.data/05.SplitGenes


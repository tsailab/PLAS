#!/bin/bash

##Additional work; check mode based on input files

#########################################################
#########derived from 01.folder.IDConverter.pl###########
#########################################################

echo 'Running 01.folder.IDConverter.pl ....' >> job.monitor.txt

## read in parameters required by the script
srcfolder="01.data/02.Fasta"
mode="paired-end"
mv fastaCombinePairedEnd.* 00.script/shell.script
rm -rf 00.script/shell.script.previous
mv 00.script/shell.script 00.script/shell.script.previous
mkdir -p 00.script/shell.script
rm -f flag*

for sub in 01.data/01.Fastq/*; do
  if [ -d "$sub" ]; then
        sub=$(basename ${sub})
        if [ $mode == "paired-end" ]; then
                time cat $srcfolder/$sub/$sub.R1.fasta_pairs_R1.fasta | awk '{if(NR%2==0){print $_} else{print ">"NR/2+0.5}}' > $srcfolder/$sub/$sub.R1.fasta_simple.fasta
                time cat $srcfolder/$sub/$sub.R2.fasta_pairs_R2.fasta | awk '{if(NR%2==0){print $_} else{print ">"NR/2+0.5}}' > $srcfolder/$sub/$sub.R2.fasta_simple.fasta
                time cat $srcfolder/$sub/$sub.R1.fasta_singles.fasta | awk '{if(NR%2==0){print $_} else{print ">"NR/2+0.5}}' > $srcfolder/$sub/$sub.singles.fasta_simple.fasta
    elif [ $mode == "single-end" ]; then
                cat $srcfolder/$sub/$sub.fasta | awk '{if(NR%2==0){print $_} else{print ">"NR/2+0.5}}' > $srcfolder/$sub/$sub.simple.fasta
        else
                echo 'Error: no read mode available!';
                exit 2
        fi
  fi
done

echo 'Finished converting folder ID!' >> job.monitor.txt
chmod 755 -R 01.data/02.Fasta
chmod 755 -R 00.script

:>01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta
:>01.data/05.SplitGenes/03.Full.Length/full.length.contigs.prot.fasta

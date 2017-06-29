#!/bin/bash

##Additional work; check mode based on input files
mode="paired-end"

#########################################################
#####derived from 040.folder.retrievebowtie.reads.pl#####
#########################################################
b=0
srcfolder="03.blast/03.bowtie.nucl/run.$b"
tgtfolder="04.retrieve.reads/03.bowtie.nucl/run.$b"
seqtype="nucl"
logfolder="bowtie.log/bowtie.run.$b"
blocksize="1000"
scale="genome"
sleeptime="20"
thread="24"

for chk in 01.data/05.SplitGenes/01.Protein/run.0; do
  if [ -d "$chk" ]; then
  chk=$(basename ${chk}) 
      for sam in $work/01.data/02.fasta; do
                if [ -d "$sam" ]; then
		sam=$(basename ${sam})
                                if [ -f "bowtie.$seqtype.$chk.$sam.sh.o*" && -f "bowtie.$seqtype.$chk.$sam.sh.e*" && -f "00.script/$logfolder/$chk.$sam.done.*" ]; then
                                        grep -E 'ERROR|Error|error' "bowtie.$seqtype.$chk.$sam.sh.e*" >> "00.script/$logfolder/summary.error.log"
                                else
                                        sleep $sleeptime
                                fi
                fi
        done
  fi
 done

        if [ -s "00.script/$logfolder/summary.error.log" ]; then
                echo "Error in retrieving bowtie reads! Check 00.script/$logfolder/summary.error.log!" >> job.monitor.txt
        	exit
	else
                echo "All jobs finished without error!" >> job.monitor.txt
        fi

        sleep "20"

        chmod 777 -R 00.script/$logfolder
        chmod 777 -R $srcfolder
        rm -rf 00.script/04.retrieve.script/run.$b
        mkdir -p 00.script/04.retrieve.script/run.$b
        rm -rf $tgtfolder
	echo "Shuffled around the files and folders correctly" >> job.monitor.txt
for sub in 03.blast/03.bowtie.nucl/run.$b/*; do
	if [ -d "$sub" ]; then
	sub=$(basename ${sub})
                mkdir -p $tgtfolder/$sub
		echo "Made directory $tgtfolder/$sub" >> job.monitor.txt
	
		for sam in 01.data/02.Fasta/*; do
			if [ -d "$sam" ]; then
			echo $sam >> job.monitor.txt
			sam=$(basename ${sam})
				if [ $mode == "paired-end" ]; then	
					echo "Starting paired-end for $sam/$sub" >> job.monitor.txt
					perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.R1.tab 01.data/02.Fasta/$sam/$sam.R1.fasta_pairs_R1.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta                        
					perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.R2.tab 01.data/02.Fasta/$sam/$sam.R2.fasta_pairs_R2.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R2.fasta                        
					perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.single.tab 01.data/02.Fasta/$sam/$sam.R1.fasta_singles.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta
					echo "Starting sed for $sam/$sub" >> job.monitor.txt
					sed -i '/^>[A-Za-z0-9_]*/s/$/\/1/' $tgtfolder/$sub/retrieved.$sub.R1.fasta
					sed -i '/^>[A-Za-z0-9_]*/s/$/\/2/' $tgtfolder/$sub/retrieved.$sub.R2.fasta
        			elif [ $mode == "single-end" ]; then                        
					perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.tab 01.data/02.Fasta/$sam/$sam.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.fasta                
				else
					echo "Error: No read type selected!" >> job.monitor.txt
				fi	
			fi
		done
		touch "00.script/04.retrieve.script/run.$b/$sub.done.log"
	fi
done
	echo 'Finished 040.folder.retrievebowtie.reads.pl!' >> job.monitor.txt
        chmod 777 -R 00.script/04.retrieve.script/run.$b


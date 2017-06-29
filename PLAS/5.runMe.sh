#!/bin/bash

##Additional work; check mode based on input files

#########################################################
###########derived from 03.diamond.folder.pl#############
#########################################################
	mode="paired-end"
	b="0"
	evalue=1e-3
        qryfolder="01.data/02.Fasta"
        dbfolder="01.data/05.SplitGenes/01.Protein/run.$b"
        tgtfolder="03.blast/03.bowtie.nucl/run.$b"
        seqtype="nucl"
        logfolder="bowtie.log/bowtie.run.$b"
        echo 'Running 03.diamond.folder.pl ....' >> job.monitor.txt
        rm -rf 00.script/$logfolder
        mkdir -p 00.script/$logfolder

	##MPI INTEGRATION HER
	for db in $dbfolder/*; do
	if [ -d "$db" ]; then
	db=$(basename ${db})
	echo "DB: $db" >> job.monitor.txt
		for sub in $qryfolder/*; do
		if [ -d "$sub" ]; then 
		sub=$(basename ${sub})
		echo "Sub: $sub" >> job.monitor.txt
		mkdir -p $tgtfolder/$db
			if [ "$mode" == "paired-end" ]; then
				echo "Running diamond blastx 1" >> job.monitor.txt
				$WORK/bin/diamond blastx -k 1 -e $evalue -d $dbfolder/$db/$db -q $qryfolder/$sub/$sub.R1.fasta_simple.fasta -a $tgtfolder/$db/bowtie.out.$db.$sub.R1
				$WORK/bin/diamond view -a $tgtfolder/$db/bowtie.out.$db.$sub.R1.daa -o $tgtfolder/$db/bowtie.out.$db.$sub.R1.tab -f tab	
				rm -f $tgtfolder/$db/bowtie.out.$db.$sub.R1.daa

				echo "Running diamond blastx 2" >> job.monitor.txt
				$WORK/bin/diamond blastx -k 1 -e $evalue -d $dbfolder/$db/$db -q $qryfolder/$sub/$sub.R2.fasta_simple.fasta -a $tgtfolder/$db/bowtie.out.$db.$sub.R2			
				$WORK/bin/diamond view -a $tgtfolder/$db/bowtie.out.$db.$sub.R2.daa -o $tgtfolder/$db/bowtie.out.$db.$sub.R2.tab -f tab				
				rm -f $tgtfolder/$db/bowtie.out.$db.$sub.R2.daa

	if [ -s "$qryfolder/$sub/$sub.singles.fasta_simple.fasta" ]; then
			$WORK/bin/diamond blastx -k 1 -e $evalue -d $dbfolder/$db/$db -q $qryfolder/$sub/$sub.singles.fasta_simple.fasta -a $tgtfolder/$db/bowtie.out.$db.$sub.single
			$WORK/bin/diamond view -a $tgtfolder/$db/bowtie.out.$db.$sub.single.daa -o $tgtfolder/$db/bowtie.out.$db.$sub.single.tab -f tab
			rm -f $tgtfolder/$db/bowtie.out.$db.$sub.single.daa
				else
					:> $tgtfolder/$db/bowtie.out.$db.$sub.single.tab
				fi
			elif [ "$mode" == "single-end" ]; then
				$WORK/bin/diamond blastx -k 1 -e $evalue -d $dbfolder/$db/$db -q $qryfolder/$sub/$sub.simple.fasta -a $tgtfolder/$db/bowtie.out.$db.$sub
				$WORK/bin/diamond view -a $tgtfolder/$db/bowtie.out.$db.$sub.daa -o $tgtfolder/$db/bowtie.out.$db.$sub.tab -f tab
				rm -f $tgtfolder/$db/bowtie.out.$db.$sub.daa
			else
				echo 'Error: No read mode available!' >> job.monitor.txt
			fi
			touch 00.script/$logfolder/$db.$sub.done.log
		fi
		done
	fi
	done		

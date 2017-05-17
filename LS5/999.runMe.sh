#!/bin/bash
#arguments here

##Additional work; check mode based on input files
mode="paired-end"
#########################################################
#######derived from 02.makeblastdb.folder.pl#############
#########################################################

echo 'Loading BLAST module...' >> job.monitor.txt
module load blast

echo 'Creating databases...' >> job.monitor.txt

##Get and loop through subdirectories
for sub in $WORK/01.data/05.SplitGenes/01.Protein/run.0/*; do	##Make sure $WORK variable works for CyVerse submission
  if [ -d "$sub" ]; then
  
	#Make blast prot db
    makeblastdb -in 01.data/05.SplitGenes/01.Protein/run.0/$sub/$sub.fasta -dbtype prot
	echo 'Created blast prot db in $sub...' >> job.monitor.txt
	
	#Make DIAMOND prot db
	$WORK/bin/diamond makedb --in 01.data/05.SplitGenes/01.Protein/run.0/$sub/$sub.fasta -d 01.data/05.SplitGenes/01.Protein/run.0/$sub/$sub
	echo 'Created diamond db in $sub...' >> job.monitor.txt

  fi
done

for sub in $WORK/01.data/05.SplitGenes/02.Transcript/run.0/*; do
  if [ -d "$sub" ]; then
	#Make blast nucl db
    makeblastdb -in 01.data/05.SplitGenes/02.Transcript/run.0/$sub/$sub.fasta -dbtype nucl
	echo 'Created blast nucl db in $sub...' >> job.monitor.txt

  fi
done

echo 'Finished creating databases!' >> job.monitor.txt
chmod 755 -R 00.script
chmod 755 -R 01.data/05.SplitGenes

#########################################################
#####derived from 01.folder.fastaCombinePairedEnd.pl#####
#########################################################

echo 'Separating paired-end and single reads...' >> job.monitor.txt
module load python

## start running the script
mv fastx.fastq2fasta.* 00.script/shell.script
rm -rf 00.script/shell.script.previous
mv 00.script/shell.script 00.script/shell.script.previous
mkdir -p 00.script/shell.script

##Get and loop through subdirectories
for sub in $WORK/01.data/01.Fastq; do
  if [ -d "$sub" ]; then
	echo 'Running paired end script for $sub' >> job.monitor.txt
	time python2.7 00.script/01.fastaCombinePairedEnd.py 01.data/01.Fastq/$sub/$sub.R1.fastq 01.data/01.Fastq/$sub/$sub.R2.fastq " "
  fi
done

echo 'Finished separating reads!' >> job.monitor.txt
chmod 755 -R 01.data/01.Fastq	##May be unnecessary
chmod 755 -R 01.data/02.Fasta
chmod 755 -R 00.script

#########################################################
#########derived from 01.fastq2fasta.folder.pl###########
#########################################################

echo 'Running 01.fastq2fasta.folder.pl ....' >> job.monitor.txt

srcfolder="01.data/01.Fastq"
tgtfolder="01.data/02.Fasta"

thread="1"

rm -rf 00.script/shell.script.previous
mv 00.script/shell.script 00.script/shell.script.previous
mkdir -p 00.script/shell.script

for sub in $WORK/01.data/01.Fastq; do
  if [ -d "$sub" ]; then
	mkdir -p 01.data/02.Fasta/$sub
	if [$mode -eq "paired-end"]; then
		awk '1 == (NR) % 4 || 2 == (NR) % 4' 01.data/01.Fastq/$sub/$sub.R1.fastq_pairs_R1.fastq | awk '{gsub(\"^@\", \">\", \$0); print \$0}' > 01.data/02.Fasta/$sub/$sub.R1.fasta_pairs_R1.fasta
		awk '1 == (NR) % 4 || 2 == (NR) % 4' 01.data/01.Fastq/$sub/$sub.R2.fastq_pairs_R2.fastq | awk '{gsub(\"^@\", \">\", \$0); print \$0}' > 01.data/02.Fasta/$sub/$sub.R2.fasta_pairs_R2.fasta
		awk '1 == (NR) % 4 || 2 == (NR) % 4' 01.data/01.Fastq/$sub/$sub.R1.fastq_singles.fastq | awk '{gsub(\"^@\", \">\", \$0); print \$0}' > 01.data/02.Fasta/$sub/$sub.R1.fasta_singles.fasta
		sed -i \"s/\\/1//g\" 01.data/02.Fasta/$sub/$sub.R1.fasta_singles.fasta
		sed -i \"s/\\/2//g\" 01.data/02.Fasta/$sub/$sub.R1.fasta_singles.fasta
	elif [ $mode -eq "single-end" ]; then
		awk '1 == (NR) % 4 || 2 == (NR) % 4' 01.data/01.Fastq/$sub/$sub.fastq | awk '{gsub(\"^@\", \">\", \$0); print \$0}' > 01.data/02.Fasta/$sub/$sub.fasta
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

#########################################################
#########derived from 01.folder.IDConverter.pl###########
#########################################################

echo 'Running 01.folder.IDConverter.pl ....' >> job.monitor.txt

## read in parameters required by the script
srcfolder="01.data/02.Fasta"



mv fastaCombinePairedEnd.* 00.script/shell.script
rm -rf 00.script/shell.script.previous
mv 00.script/shell.script 00.script/shell.script.previous
mkdir -p 00.script/shell.script
rm -f flag*


for sub in $WORK/01.data/01.Fastq; do
  if [ -d "$sub" ]; then
	if [ $mode -eq "paired-end" ]; then
		time cat $srcfolder/$sub/$sub.R1.fasta_pairs_R1.fasta | awk '{if(NR%2==0){print \$_} else{print \">\"NR/2+0.5}}' > $srcfolder/$sub/$sub.R1.fasta_simple.fasta
		time cat $srcfolder/$sub/$sub.R2.fasta_pairs_R2.fasta | awk '{if(NR%2==0){print \$_} else{print \">\"NR/2+0.5}}' > $srcfolder/$sub/$sub.R2.fasta_simple.fasta
		time cat $srcfolder/$sub/$sub.R1.fasta_singles.fasta | awk '{if(NR%2==0){print \$_} else{print \">\"NR/2+0.5}}' > $srcfolder/$sub/$sub.singles.fasta_simple.fasta
    elif [ $mode -eq "single-end" ]; then
		cat $srcfolder/$sub/$sub.fasta | awk '{if(NR%2==0){print \$_} else{print \">\"NR/2+0.5}}' > $srcfolder/$sub/$sub.simple.fasta
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

#########################################################
##########TURBOTURBOTURBOTURBOTURBOTURBOTURBO############
#########################################################

b=0
evalue=1e-3
while [ $b -le 14 ] #number of runs; make user configurable.
do
    echo "Run number: $b" >> job.monitor.txt

    if [ $b -eq 0 ];then

#########################################################
###########derived from 03.diamond.folder.pl#############
#########################################################

	qryfolder="01.data/02.Fasta"
	dbfolder="01.data/05.SplitGenes/01.Protein/run.$b"
	tgtfolder="03.blast/03.bowtie.nucl/run.$b"
	seqtype="nucl"
	logfolder="bowtie.log/bowtie.run.$b"
	thread="4"
	echo 'Running 03.diamond.folder.pl ....' >> job.monitor.txt			
	rm -rf 00.script/$logfolder
	mkdir -p 00.script/$logfolder	
		
		##MPI INTEGRATION HERE
		for db in $WORK/$dbfolder; do	
			if [ -d "$db" ]; then
				for sub in $WORK/$qryfolder; do
					if [ -d "$sub" ]; then
						mkdir -p $tgtfolder/$db
							if [ $mode -eq "paired-end" ]; then
								diamond blastx -p $thread -k 1 -e $evalue -d $dbfolder/$db/$db -q $qryfolder/$sub/$sub.R1.fasta_simple.fasta -o $tgtfolder/$db/bowtie.out.$db.$sub.R1.tab
								diamond blastx -p $thread -k 1 -e $evalue -d $dbfolder/$db/$db -q $qryfolder/$sub/$sub.R2.fasta_simple.fasta -o $tgtfolder/$db/bowtie.out.$db.$sub.R2.tab
								if [ -f "$qryfolder/$sub/$sub.singles.fasta_simple.fasta" ] && [ ! -s "$qryfolder/$sub/$sub.singles.fasta_simple.fasta" ]; then
									$WORK/bin/diamond blastx -p $thread -k 1 -e $evalue -d $dbfolder/$db/$db -q $qryfolder/$sub/$sub.singles.fasta_simple.fasta -o $tgtfolder/$db/bowtie.out.$db.$sub.single.tab
								else
									:> $tgtfolder/$db/bowtie.out.$db.$sub.single.tab
								fi
							elif [ $mode -eq "single-end" ]; then
								diamond blastx -p $thread -k 1 -e $evalue -d $dbfolder/$db/$db -q $qryfolder/$sub/$sub.simple.fasta -o $tgtfolder/$db/bowtie.out.$db.$sub.tab
							else
								echo 'Error: No read mode available!' >> job.monitor.txt
							fi
						touch 00.script/$logfolder/$db.$sub.done.log
					fi
				done
			fi
		done
		
#########################################################
#####derived from 040.folder.retrievebowtie.reads.pl#####
#########################################################

srcfolder = "03.blast/03.bowtie.nucl/run.$b"
tgtfolder = "04.retrieve.reads/03.bowtie.nucl/run.$b"
seqtype = "nucl"
logfolder = "bowtie.log/bowtie.run.$b"
blocksize = "1000"
scale = "genome"
platform = lc(shift @ARGV)
sleeptime = "20"
##replace $run with $b
thread = 1
	
@chks = 01.data/05.SplitGenes/01.Protein/run.0 subdirectories
@sams = 01.data/02.fasta subdirectories
@subs = $srcfolder subdirectories


for chk in $WORK/01.data/05.SplitGenes/01.Protein/run.0; do
  if [ -d "$chk" ]; then
	for sam in $work/01.data/02.fasta; do
		if [ -d "$sam" ]; then
				if [ -f "bowtie.$seqtype.$chk.$sam.sh.o*" && -f "bowtie.$seqtype.$chk.$sam.sh.e*" && -f "00.script/$logfolder/$chk.$sam.done.*" ]; then
					grep -E 'ERROR|Error|error' "bowtie.$seqtype.$chk.$sam.sh.e*" >> "00.script/$logfolder/summary.error.log"
				else
					sleep $sleeptime
				fi
		fi	
	done
  fi
 done

	if [ !-s "00.script/$logfolder/summary.error.log" ]; then
		echo "All jobs finished without error!" >> job.monitor.txt
	else
		echo "Error in retrieving bowtie reads! Check 00.script/$logfolder/summary.error.log!" >> job.monitor.txt
		exit
	fi
	sleep "20"
	
	mv bowtie.* 00.script/$logfolder
	chmod 777 -R 00.script/$logfolder
	chmod 777 -R $srcfolder
	rm -rf 00.script/04.retrieve.script/run.$b
	mkdir -p 00.script/04.retrieve.script/run.$b
	rm -rf $tgtfolder
	
	for sub in $work/03.blast/03.bowtie.nucl/run.$b; do
		mkdir -p $tgtfolder/03.blast/03.bowtie.nucl/run.$b	

		if [ $mode -eq "paired-end" ]; then	#Bash script calling perl script that makes a bash script that calls a perl script that makes a... 
			perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.R1.tab 01.data/02.Fasta/$sam/$sam.R1.fasta_pairs_R1.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta
			perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.R2.tab 01.data/02.Fasta/$sam/$sam.R2.fasta_pairs_R2.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R2.fasta
			perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.single.tab 01.data/02.Fasta/$sam/$sam.R1.fasta_singles.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta
			sed -i '/^>[A-Za-z0-9_]*/s/\$/\\/1/' $tgtfolder/$sub/retrieved.$sub.R1.fasta
			sed -i '/^>[A-Za-z0-9_]*/s/\$/\\/2/' $tgtfolder/$sub/retrieved.$sub.R2.fasta
		elif [ $mode -eq "single-end" ]; then
			perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.tab 01.data/02.Fasta/$sam/$sam.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.fasta
		fi
		
	touch 00.script/04.retrieve.script/run.$b/$sub.done.log	
	echo 'Finished 040.folder.retrievebowtie.reads.pl!' >> job.monitor.txt
	chmod 777 -R 00.script/04.retrieve.script/run.$b
	
#########################################################
######TURBOTURBOTURBOTURBOTURBOTURBOTURBOTURBOTURBO######
#########################################################

	else

#########################################################	
##############021.makebowtiedb.folder.pl#################
#########################################################
	module load bowtie2
	
	tgtfolder="01.data/05.SplitGenes/02.Transcript/run.$b"
	sleeptime="10"
	let prerun="$b-1"
	thread="1"
	
	while [ "$b" -ge 0 ]; then
		if [ ! -f $WORK/00.script/10.transfer.script/run.$prerun/transfer.saturate.seq.log ]; then
			sleep $sleeptime
		elif [ -s 00.script/10.transfer.script/run.$prerun/summary.error.log ]; then
			echo "Error: Something went wrong. See 00.script/10.transfer.script/run.$prerun/summary.error.log"
		else
			echo "All jobs have finished successfully!" >> job.monitor.txt
		fi
	done
	
	rm -rf 00.script/02.makebowtiedb.script/run.$b
	mkdir -p 00.script/02.makebowtiedb.script/run.$b
	chmod 755 -R $tgtfolder
	
	for sub in $tgtfolder; do
		if [ -d "$sub" ]; then
			cd $tgtfolder/$sub
			time bowtie2-build -f -q $sub.fasta $sub
			cd ../../../../../
			touch 00.script/02.makebowtiedb.script/run.$b/$sub.done.log
		fi
	done
	
echo 'Finished 021.makeblastdb.folder.pl!' >> job.monitor.txt				
			
#########################################################	
##################03.bowtie.folder.pl####################
#########################################################			
			
	qryfolder="01.data/02.Fasta"
	dbfolder="01.data/05.SplitGenes/02.Transcript/run.$b"
	tgtfolder="03.blast/03.bowtie.nucl/run.$b"
	reffolder="01.data/05.SplitGenes/01.Protein/run.0"
	seqtype="nucl"
	logfolder="bowtie.log/bowtie.run.$b"
	sleeptime="10"
	run=$b
	thread="4"

	for chk in $WORK/$reffolder; do
		if [ -d "$chk" ]; then
			if [ -f "makebowtiedb.$chk.sh.o*" && -f "makebowtiedb.$chk.sh.e*" && -f "00.script/02.makebowtiedb.script/run.$b/$chk.done.*" ]; then
				grep -E 'ERROR|Error|error' "makebowtiedb.$chk.sh.e*" >> "00.script/$logfolder/summary.error.log"
			fi	
		fi
	done
	
	if [ !-s "00.script/$logfolder/summary.error.log" ]; then
		echo "All jobs finished without error!" >> job.monitor.txt
		for chk in $WORK/$reffolder; do
		if [ -d "$chk" ]; then
			if [ !-s "$dbfolder/$chk/$chk.1.bt2" || !-s "$dbfolder/$chk/$chk.rev.1.bt2" ]; then
				echo "No output file for $chk!" >> job.monitor.txt
			fi
		fi
		done
	else
		echo "Error in making bowtie databases! Check 00.script/$logfolder/summary.error.log!" >> job.monitor.txt
		exit
	fi
	sleep "20"
	
	mv makebowtiedb.* 00.script/02.makebowtiedb.script/run.$b
	chmod 777 -R 00.script/02.makebowtiedb.script/run.$b
	chmod 777 -R $dbfolder
	rm -rf 00.script/$logfolder
	mkdir -p 00.script/$logfolder
	
	for db in "$WORK/$dbfolder"; do
	if [ -d "$db" ]; then
		if [ $mode -eq "Paired-end" ]; then
			time bowtie2 -f -x $dbfolder/$db/$db $unmap -p $thread --local -1 "$qryfolder/$sub/$sub.R1.fasta_pairs_R1.fasta" -2 "$qryfolder/$sub/$sub.R2.fasta_pairs_R2.fasta" -U "$qryfolder/$sub/$sub.R1.fasta_singles.fasta" -S $tgtfolder/$db/bowtie.out.$db.sam
		elif [ $mode -q "Single-end" ]; then
			time bowtie2 -f -x $dbfolder/$db/$db $unmap -p $thread --local -U "$qryfolder/$sub/$sub.R1.fasta_pairs_R1.fasta" -S $tgtfolder/$db/bowtie.out.$db.sam
		else
			echo "No read mode selected!" >> job.monitor.txt
			exit
		fi
	
	touch "00.script/$logfolder/$db.done.log"
	fi
	done
	
	echo 'Finished 03.bowtie.folder.pl!' >> job.monitor.txt
	chmod 777 -R 00.script/$logfolder
	
#########################################################	
###########04.folder.retrievebowtie.reads.pl#############	##Remember to replace "run" with "b"
#########################################################		

srcfolder="03.blast/03.bowtie.nucl/run.$b"
tgtfolder="04.retrieve.reads/03.bowtie.nucl/run.$b"
seqtype="nucl"
logfolder="bowtie.log/bowtie.run.$b"
sleeptime="20"
unmap="map"
if [ $b -eq "0" ]; then
	unmap="both"
fi
#ERROR CHECKING

mv bowtie.* 00.script/$logfolder
chmod 777 -R 00.script/$logfolder
chmod 777 -R $srcfolder

rm -rf 00.script/04.retrieve.script/run.$run
mkdir -p 00.script/04.retrieve.script/run.$run

for sub in $WORK/$srcfolder; do
	if [ -d "$sub" ]; then
		mkdir -p $tgtfolder/$sub\n
		
		if [ $mode -eq "paired-end" ]; then
			perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.sam $unmap $tgtfolder/$sub/retrieved.$sub.R1.fasta $tgtfolder/$sub/unmap.$sub.R1.fasta $tgtfolder/$sub/retrieved.$sub.R2.fasta $tgtfolder/$sub/unmap.$sub.R2.fasta
		elif [ $mode -eq "single-end" ]; then
			perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.sam $unmap $tgtfolder/$sub/retrieved.$sub.fasta $tgtfolder/$sub/unmap.$sub.R1.fasta
		else	
			echo "Error, specify read mode" >> job.monitor.text
			exit
		fi
		
		if [ -e "$srcfolder/$sub/bowtie.out.$sub.long.sam" && -s "$srcfolder/$sub/bowtie.out.$sub.long.sam" ]; then	##Check this comparison
			perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.long.sam $unmap $tgtfolder/$sub/retrieved.$sub.long.fasta $tgtfolder/$sub/unmap.$sub.long.fasta
		fi
		touch 00.script/04.retrieve.script/run.$run/$sub.done.log
	fi
done

echo 'Finished 04.folder.retrievebowtie.reads.pl!' >> job.monitor.txt
chmod 777 -R 00.script/04.retrieve.script/run.$run

	
#########################################################	
###########END OF TURBOTURBOTURBOTURBOTURBO$#############
#########################################################	

	fi
	
#########################################################	
#############06.assembly.trinity.folder.pl###############
#########################################################
	
04.retrieve.reads/03.bowtie.nucl/run.$b 
06.assembly/03.bowtie.nucl/run.$b 
genome 
$mode 
$platform 
20

#########################################################	
##############06.truncate.header.folder.pl###############
#########################################################
	
06.assembly/03.bowtie.nucl/run.$b 
$platform 
20

	
#########################################################	
###################07.blastx.back.pl#####################
#########################################################
echo 'Running 07.blastx.back.pl ....' >> job.monitor.txt
module load blast
srcfolder="06.assembly/03.bowtie.nucl/run.$b"
dbfolder="01.data/05.SplitGenes/01.Protein/run.0"
tgtfolder="07.map.back/03.bowtie.nucl/run.$b"
reffolder="01.data/05.SplitGenes/01.Protein/run.0"
thread="1"

###ERROR CHECKING

mv truncate.header.* 00.script/06.truncate.script/run.$run
chmod 777 -R 00.script/06.truncate.script/run.$run
chmod 777 -R $srcfolder
rm -rf 00.script/07.blastx.script/run.$run
mkdir -p 00.script/07.blastx.script/run.$run

for sub in "$WORK/$srcfolder"; do
	if [ -d "$sub" ]; then
		mkdir -p $tgtfolder/$sub
		blastx -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.out -evalue 1e-5  -outfmt 6 -num_threads 1 -max_target_seqs 1
		blastx -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.xml.out -evalue 1e-5  -outfmt 5 -num_threads 1 -max_target_seqs 1
		touch 00.script/07.blastx.script/run.$run/$sub.done.log
		
		
		
	fi
done

echo 'Finished 07.blastx.back.pl!' >> job.monitor.txt
chmod 777 -R 00.script/07.blastx.script/run.$run

#########################################################	
###################07.blastn.back.pl#####################
#########################################################

echo 'Running 07.blastx.back.pl ....' >> job.monitor.txt

srcfolder="06.assembly/03.bowtie.nucl/run.$b"
dbfolder="01.data/05.SplitGenes/02.Transcript/run.0"
tgtfolder="07.map.back/02.blastn/run.$b"
thread="1"

###ERROR CHECKING

rm -rf 00.script/07.blastn.script/run.$run
mkdir -p 00.script/07.blastn.script/run.$run

for sub in "$WORK/$srcfolder"; do
	if [ -d "$sub" ]; then
		blastn -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.out -evalue 1e-10  -outfmt 6 -num_threads 1 -max_target_seqs 1
		blastn -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.xml.out -evalue 1e-10  -outfmt 5 -num_threads 1 -max_target_seqs 1
	fi
done

echo 'Finished 07.blastx.back.pl!' >> job.monitor.txt
chmod 777 -R 00.script/07.blastn.script/run.$run

    c=`expr $b + 1`

    time perl 00.script/100.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.$b 07.map.back/03.bowtie.nucl/run.$b 01.data/05.SplitGenes/02.Transcript/run.$c 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 
    time perl 00.script/10.folder.detect.full.length.seq.pl 07.map.back/03.bowtie.nucl/run.$b 06.assembly/03.bowtie.nucl/run.$b 01.data/05.SplitGenes/01.Protein/run.0 01.data/05.SplitGenes/03.Full.Length/run.$c pct 0.1
    cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.nucl.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta
    cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.prot.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.prot.fasta

    a=`expr $a + 1`	
	
done

#### summarize all runs
grep ">" 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta > 01.data/05.SplitGenes/03.Full.Length/count1
sed -i "s/>//" 01.data/05.SplitGenes/03.Full.Length/count1
perl 00.script/b3.full.length.format.pl 01.data/04.GeneOfInterest/GeneID.v1.txt 01.data/05.SplitGenes/03.Full.Length/count1 01.data/05.SplitGenes/03.Full.Length/count2 01.data/05.SplitGenes/03.Full.Length/count3
wait

#### assemble unmapped reads
cp 01.data/05.SplitGenes/03.Full.Length/count3 08.full.length/
cp 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta 08.full.length/
time perl 00.script/c9.get.full.length.seq.pl 08.full.length/count3 08.full.length/full.length.contigs.nucl.fasta 08.full.length/Final.v1.fasta

blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 32 -max_target_seqs 1

makeblastdb -in 08.full.length/Final.v1.fasta -dbtype nucl

blastn -db 08.full.length/Final.v1.fasta -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5

time perl 00.script/c11.remove.redundancy.pl 08.full.length/Final.v1.blastn.xml.out 08.full.length/Final.v1.fasta 08.full.length/Final.v2.fasta self 08.full.length/Final.v1.ptr.blastx.out > 08.full.length/record.run1
wait

makeblastdb -in 08.full.length/Final.v2.fasta -dbtype nucl

blastn -db 08.full.length/Final.v2.fasta -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5

blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 32 -max_target_seqs 1

time perl 00.script/101.transfer.saturate.seq.pl 08.full.length/Final.v2.ptr.blastx.out 08.full.length/Final.v2.fasta 01.data/00.PriorData/ptr.proteome.fa 08.full.length Final.v2.ptr pct 0.02
wait

cd 08.full.length/
ln -sf Final.v2.ptr.full.contigs.nucl.fasta Final.fasta
cd ../



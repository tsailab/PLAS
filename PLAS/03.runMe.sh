#PBS -S /bin/bash
#PBS -N runMe
#PBS -q highmem_q
#PBS -l nodes=1:ppn=8
#PBS -l mem=60gb
#PBS -l walltime=12:00:00
#PBS -d ./

cd $PBS_O_WORKDIR
sub=${SUB}
mode="paired-end"
logfolder="00.script/01.log"
run="0"

error_check() {
	if [[ $? -ne 0 ]]; then
	echo "$1" >> 00.script/01.log/job.monitor_$sub.txt
	exit $?
	fi
return 0
}

# Load required modules
module load BLAST+/2.6.0-foss-2016b-Python-2.7.14
module load Python/2.7.14-foss-2016b
module load diamond/1.0
module load Bowtie2/2.3.3-foss-2016b
module load Trinity/2.6.6-foss-2016b

# Create BLAST databases for protein and transcripts 


echo "Creating database for $sub..." >> $logfolder/job.monitor_$sub.txt
makeblastdb -in 01.data/05.SplitGenes/01.Protein/run.0/$sub/$sub.fasta -dbtype prot
error_check "Failed to make blast prot db for $sub, check runMe.sh line "

diamond makedb --in 01.data/05.SplitGenes/01.Protein/run.0/$sub/$sub.fasta -d 01.data/05.SplitGenes/01.Protein/run.0/$sub/$sub
error_check "Failed to make diamond prot db for $sub, check runMe.sh line "

makeblastdb -in 01.data/05.SplitGenes/02.Transcript/run.0/$sub/$sub.fasta -dbtype nucl

echo 'Finished creating databases!' >> $logfolder/job.monitor_$sub.txt


# Divided all proteins or transcripts of reference genome into different meta-groups by DIAMOND 

mode="paired-end"
b="0"
evalue=1e-3
qryfolder="01.data/02.Fasta"
subfolder="01.data/05.SplitGenes/01.Protein/run.$run"
tgtfolder="03.blast/03.bowtie.nucl/run.$run"
seqtype="nucl"
echo "DB: $sub" >> $logfolder/job.monitor_$sub.txt
	for qry in $qryfolder/*; do
	if [ -d "$qry" ]; then 
	qry=$(basename ${qry})
	echo "qry: $qry" >> $logfolder/job.monitor_$sub.txt
	mkdir -p $tgtfolder/$sub
	
		if [ "$mode" == "paired-end" ]; then
			echo "Running diamond blastx 1" >> $logfolder/job.monitor_$sub.txt
			diamond blastx -k 1 -e $evalue -d $subfolder/$sub/$sub -q $qryfolder/$qry/$qry.R1.fasta_simple.fasta -a $tgtfolder/$sub/bowtie.out.$sub.$qry.R1
			diamond view -a $tgtfolder/$sub/bowtie.out.$sub.$qry.R1.daa -o $tgtfolder/$sub/bowtie.out.$sub.$qry.R1.tab -f tab	
			rm -f $tgtfolder/$sub/bowtie.out.$sub.$qry.R1.daa

			echo "Running diamond blastx 2" >> $logfolder/job.monitor_$qry.txt
			diamond blastx -k 1 -e $evalue -d $subfolder/$sub/$sub -q $qryfolder/$qry/$qry.R2.fasta_simple.fasta -a $tgtfolder/$sub/bowtie.out.$sub.$qry.R2			
			diamond view -a $tgtfolder/$sub/bowtie.out.$sub.$qry.R2.daa -o $tgtfolder/$sub/bowtie.out.$sub.$qry.R2.tab -f tab				
			rm -f $tgtfolder/$sub/bowtie.out.$sub.$qry.R2.daa

			if [ -s "$qryfolder/$qry/$qry.singles.fasta_simple.fasta" ]; then
				diamond blastx -k 1 -e $evalue -d $subfolder/$sub/$sub -q $qryfolder/$qry/$qry.singles.fasta_simple.fasta -a $tgtfolder/$sub/bowtie.out.$sub.$qry.single
				diamond view -a $tgtfolder/$sub/bowtie.out.$sub.$qry.single.daa -o $tgtfolder/$sub/bowtie.out.$sub.$qry.single.tab -f tab
				rm -f $tgtfolder/$sub/bowtie.out.$sub.$qry.single.daa
			else
				:> $tgtfolder/$sub/bowtie.out.$sub.$qry.single.tab
			fi
		
		elif [ "$mode" == "single-end" ]; then
			diamond blastx -k 1 -e $evalue -d $subfolder/$sub/$sub -q $qryfolder/$qry/$qry.simple.fasta -a $tgtfolder/$sub/bowtie.out.$sub.$qry
			diamond view -a $tgtfolder/$sub/bowtie.out.$sub.$qry.daa -o $tgtfolder/$sub/bowtie.out.$sub.$qry.tab -f tab
			rm -f $tgtfolder/$sub/bowtie.out.$sub.$qry.daa
		else
			echo 'Error: No read mode available!' >> $logfolder/job.monitor_$sub.txt
		fi
		
	fi
	done	

# Retrieve reads from bowtie results
srcfolder="03.blast/03.bowtie.nucl/run.$run"
tgtfolder="04.retrieve.reads/03.bowtie.nucl/run.$run"
seqtype="nucl"
blocksize="1000"
scale="genome"
sleeptime="20"
thread="24"


chmod 777 -R $logfolder
mkdir -p $tgtfolder/$sub
	
for sam in 01.data/02.Fasta/*; do
	if [ -d "$sam" ]; then
	sam=$(basename ${sam})
		if [ $mode == "paired-end" ]; then	
			echo "Starting paired-end for $sam/$sub" >> $logfolder/job.monitor_$sub.txt
			perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.R1.tab 01.data/02.Fasta/$sam/$sam.R1.fasta_pairs_R1.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta                        
			perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.R2.tab 01.data/02.Fasta/$sam/$sam.R2.fasta_pairs_R2.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R2.fasta                        
			perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.single.tab 01.data/02.Fasta/$sam/$sam.R1.fasta_singles.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta
			sed -i '/^>[A-Za-z0-9_]*/s/$/\/1/' $tgtfolder/$sub/retrieved.$sub.R1.fasta
			sed -i '/^>[A-Za-z0-9_]*/s/$/\/2/' $tgtfolder/$sub/retrieved.$sub.R2.fasta
			elif [ $mode == "single-end" ]; then                        
			perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.tab 01.data/02.Fasta/$sam/$sam.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.fasta                
		else
			echo "Error: No read type selected!" >> $logfolder/job.monitor_$sub.txt
		fi	
	fi
done

## Script ends

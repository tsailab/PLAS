#!/bin/bash
#SBATCH -J PLAS_runMe7
#SBATCH -o PLAS_runMe7.o%j
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal
#SBATCH -t 04:00:00

#SBATCH --mail-user=sjq28742@uga.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
module load bioperl
module load blast
module load trinityrnaseq
##

mode="paired-end"
b="3"
echo "Initial value: $b" >> checkit.txt
evalue=1e-3
touch job.monitor.txt
chmod 755 job.monitor.txt

error_check() {
        if [[ $? -ne 0 ]]; then
        echo "$1" >> job.monitor.txt
        exit $?
        fi
return 0
}

while [ $b -le 3 ]; do #number of runs; make user configurable.
   
    echo "Run number: $b" >> job.monitor.txt
if [ $b -gt "0" ]; then
    tgtfolder="01.data/05.SplitGenes/02.Transcript/run.$b"
    sleeptime="10"
    let prerun="$b-1"
    thread="24"
    module load bowtie/2.2.6
    mkdir $tgtfolder
    chmod 755 -R $tgtfolder        
        
	for sub in $tgtfolder/*; do
                if [ -d "$sub" ]; then
		sub=$(basename ${sub})
                	cd $tgtfolder/$sub
			bowtie2-build -f -q $sub.fasta $sub                        
			cd ../../../../../
		fi
        done

	echo "Finished 021.makebowtiedb.folder.pl!" >> job.monitor.txt

################
    qryfolder="01.data/02.Fasta"
    dbfolder="01.data/05.SplitGenes/02.Transcript/run.$b"
    tgtfolder="03.blast/03.bowtie.nucl/run.$b"
    reffolder="01.data/05.SplitGenes/01.Protein/run.0"
    seqtype="nucl"
    logfolder="bowtie.log/bowtie.run.$b"
    thread="24"
    R1=()
    R2=()
    R3=()
    str2=""
    str3=""

   	mkdir $tgtfolder 
	for chk in $reffolder/*; do
                if [ -d "$chk" ]; then
		chk=$(basename ${chk})
                        if [[ ! -s "$dbfolder/$chk/$chk.1.bt2" || ! -s "$dbfolder/$chk/$chk.rev.1.bt2" ]]; then
                                echo "No output file for $chk!" >> job.monitor.txt
			fi
                fi
    done

    #chmod 777 -R $dbfolder

    for db in $dbfolder/*; do
        if [ -d "$db" ]; then
	    db=$(basename ${db})
		mkdir -p $tgtfolder/$db

		if [ "$mode" == "paired-end" ]; then
		##Construct R1 input string
			str1=""
			R1=()
			for sub in $qryfolder/*; do cd $sub
				for item in ./*"R1"."fasta"; do [ -f "$item" ] || continue
				item=$(basename ${item}); R1+=("$sub/$item")
				done; cd ../../../
			done

			for elem in "${R1[@]}"; do str1="$str1,$elem"; 
			done
			str1=${str1#?}
			
		##Construct R2 input string
                        str2=""
			R2=()
			for sub in $qryfolder/*; do cd $sub
                                for item in ./*"R2"."fasta"; do [ -f "$item" ] || continue
                                item=$(basename ${item}); R2+=("$sub/$item")
                                done; cd ../../../
                        done

                        for elem in "${R2[@]}"; do str2="$str2,$elem";
                        done
                        str2=${str2#?}
		
		##Construct singles input string
                        str3=""
			R3=()
			for sub in $qryfolder/*; do cd $sub
                                for item in ./*"singles"."fasta"; do [ -f "$item" ] || continue
                                item=$(basename ${item}); R3+=("$sub/$item")
                                done; cd ../../../
                        done

                        for elem in "${R3[@]}"; do str3="$str3,$elem";
                        done
                        str3=${str3#?}
							#-unmap
			time bowtie2 -f -x $dbfolder/$db/$db -p $thread --local -1 "$str1" -2 "$str2" -S $tgtfolder/$db/bowtie.out.$db.sam
       

	         elif [ "$mode" == "Single-end" ]; then
		##Construct R1 input string
                        str1=""
			R1=()
			for sub in $qryfolder/*; do cd $sub
                                for item in ./*"R1"."fasta"; do [ -f "$item" ] || continue
                                item=$(basename ${item}); R1+=("$sub/$item")
                                done; cd ../../../
                        done

                        for elem in "${R1[@]}"; do str1="$str1,$elem";
                        done
                        str1=${str1#?}
							#$unmap
                        time bowtie2 -f -x $dbfolder/$db/$db -p $thread --local -U "$str1" -S $tgtfolder/$db/bowtie.out.$db.sam
                else
                        echo "No read mode selected!" >> job.monitor.txt
                        exit
                fi
	fi
    done

    echo 'Finished 03.bowtie.folder.pl!' >> job.monitor.txt
						
#############################

	srcfolder="03.blast/03.bowtie.nucl/run.$b"
	tgtfolder="04.retrieve.reads/03.bowtie.nucl/run.$b"
	seqtype="nucl"
	logfolder="bowtie.log/bowtie.run.$b"
	sleeptime="20"
	unmap="map"
	if [ $b -eq "0" ]; then
		echo "Somethin screwy going on here" >> job.monitor.txt
		exit
	fi
	mkdir $tgtfolder	
	chmod 777 -R $srcfolder

	for sub in $srcfolder/*; do
		if [ -d "$sub" ]; then
		sub=$(basename ${sub})
            mkdir -p $tgtfolder/$sub
			if [ "$mode" == "paired-end" ]; then
				perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.sam $unmap $tgtfolder/$sub/retrieved.$sub.R1.fasta $tgtfolder/$sub/unmap.$sub.R1.fasta $tgtfolder/$sub/retrieved.$sub.R2.fasta $tgtfolder/$sub/unmap.$sub.R2.fasta
				error_check "Bowtie broke at $sub"
			elif [ "$mode" == "single-end" ]; then                        
				perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.sam $unmap $tgtfolder/$sub/retrieved.$sub.fasta $tgtfolder/$sub/unmap.$sub.R1.fasta
			else
				echo "Error, specify read mode" >> job.monitor.txt
				exit
            fi
			touch "00.script/04.retrieve.script/run.$run/$sub.done.log"
        fi
	done

echo 'Finished 04.folder.retrievebowtie.reads.pl!' >> job.monitor.txt
chmod 777 -R 00.script/04.retrieve.script/run.$run

fi

#######################################
echo "running trinity.runMe.sh ..." >> job.monitor.txt

srcfolder="04.retrieve.reads/03.bowtie.nucl/run.$b"

chmod 777 -R $srcfolder
chmod 777 trinity.runMe.sh
for sub in $srcfolder/*; do
	if [ -d "$sub" ]; then
	sub=$(basename ${sub})
		./trinity.runMe.sh "$b" "$sub"
	fi
done
		echo "Finished running trinity.runMe.sh!" >> job.monitor.txt
#######################################
###06.truncate.header.folder.pl
srcfolder="06.assembly/03.bowtie.nucl/run.$b"
sleeptime="20"
run=$b
thread="24"

echo "Running 06.truncate.header.folder.pl ..." >> job.monitor.txt

for sub in $srcfolder/*; do
if [ -d $sub ]; then
sub=$(basename ${sub})
	time perl 00.script/06.truncate.header.pl $srcfolder/$sub/Trinity.fasta $srcfolder/$sub/Trinity.new.fasta
fi
done
echo "Finished running 06.truncate.header.folder.pl!" >> job.monitor.txt
#######################################

echo "Running 07.blastx.back.pl ..." >> job.monitor.txt
srcfolder="06.assembly/03.bowtie.nucl/run.$b"
dbfolder="01.data/05.SplitGenes/01.Protein/run.0"
tgtfolder="07.map.back/03.bowtie.nucl/run.$b"
reffolder="01.data/05.SplitGenes/01.Protein/run.0"
thread="24"

mkdir $tgtfolder

for sub in $srcfolder/*; do
	if [ -d "$sub" ]; then
	sub=$(basename ${sub})
	newsub=${sub%????????}
		mkdir -p $tgtfolder/$newsub
		blastx -db $dbfolder/$newsub/$newsub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$newsub/$newsub.contigs.blast.out -evalue 1e-5  -outfmt 6 -num_threads 1 -max_target_seqs 1
		blastx -db $dbfolder/$newsub/$newsub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$newsub/$newsub.contigs.blast.xml.out -evalue 1e-5  -outfmt 5 -num_threads 1 -max_target_seqs 1
		
	fi
done

echo 'Finished 07.blastx.back.pl!' >> job.monitor.txt

#######################################


echo "Running 07.blastn.back.pl ..." >> job.monitor.txt

srcfolder="06.assembly/03.bowtie.nucl/run.$b"
dbfolder="01.data/05.SplitGenes/02.Transcript/run.0"
tgtfolder="07.map.back/02.blastn/run.$b"
thread="24"

mkdir $tgtfolder

for sub in $srcfolder/*; do
	if [ -d "$sub" ]; then
	sub=$(basename ${sub})
	newsub=${sub%????????}
		mkdir -p $tgtfolder/$newsub
		blastn -db $dbfolder/$newsub/$newsub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$newsub/$newsub.contigs.blast.out -evalue 1e-10  -outfmt 6 -num_threads 1 -max_target_seqs 1
		blastn -db $dbfolder/$newsub/$newsub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$newsub/$newsub.contigs.blast.xml.out -evalue 1e-10  -outfmt 5 -num_threads 1 -max_target_seqs 1
	fi
done

echo "Finished 07.blastn.back.pl!" >> job.monitor.txt

#######################################

c=`expr $b + 1`

#######################################
	time perl 00.script/100.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.$b 07.map.back/03.bowtie.nucl/run.$b 01.data/05.SplitGenes/02.Transcript/run.$c 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 
#######################################
	time perl 00.script/10.folder.detect.full.length.seq.pl 07.map.back/03.bowtie.nucl/run.$b 06.assembly/03.bowtie.nucl/run.$b 01.data/05.SplitGenes/01.Protein/run.0 01.data/05.SplitGenes/03.Full.Length/run.$c pct 0.1
#######################################

cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.nucl.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta
cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.prot.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.prot.fasta
echo "Before increment: $b" >> "checkit.txt"
a=`expr $a + 1`
b=`expr $b + 1`
echo "After increment: $b" >> "checkit.txt"
done


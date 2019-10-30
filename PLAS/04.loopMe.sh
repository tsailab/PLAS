#PBS -S /bin/bash
#PBS -N loopMe
#PBS -q highmem_q
#PBS -l nodes=1:ppn=12
#PBS -l mem=48gb
#PBS -l walltime=10:00:00
#PBS -d ./

###################
###LOOP OF DEATH###
###################
cd $PBS_O_WORKDIR
run=${RUN}
sub=${SUB}

#Load required modules
module load Bowtie2/2.3.3-foss-2016b
module load BLAST+/2.6.0-foss-2016b-Python-2.7.14
module load Perl/5.20.3-foss-2016b
module load Trinity/2.6.6-foss-2016b

evalue=1e-3

error_check() {
        if [[ $? -ne 0 ]]; then
        echo $1
        exit $?
        fi
return 0
}


mode="paired-end"
if [ $run -gt "0" ]; then
	tgtfolder="01.data/05.SplitGenes/02.Transcript/run.$run"
	sleeptime="10"
	let prerun="$run-1"
    cd $tgtfolder/$sub	
    bowtie2-build -f -q $sub.fasta $sub                        
	error_check "Bowtie build failed at line 39!"
	cd ../../../../../
		
################
    qryfolder="01.data/02.Fasta"
    dbfolder="01.data/05.SplitGenes/02.Transcript/run.$run"
    tgtfolder="03.blast/03.bowtie.nucl/run.$run"
    reffolder="01.data/05.SplitGenes/01.Protein/run.0"
    seqtype="nucl"
    logfolder="bowtie.log/bowtie.run.$run"
    thread="24"
    R1=()
    R2=()
    R3=()
    str2=""
    str3=""
		
	mkdir -p $tgtfolder/$sub

		if [ "$mode" == "paired-end" ]; then
		##Construct R1 input string
			str1=""
			R1=()
			for sam in $qryfolder/*; do cd $sam
				for item in ./*"R1"."fasta"; do [ -f "$item" ] || continue
				item=$(basename ${item}); R1+=("$sam/$item")
				done; cd ../../../
			done

			for elem in "${R1[@]}"; do str1="$str1,$elem"; done
			str1=${str1#?}
			
		##Construct R2 input string
                        str2=""
			R2=()
			for sam in $qryfolder/*; do cd $sam
                                for item in ./*"R2"."fasta"; do [ -f "$item" ] || continue
                                item=$(basename ${item}); R2+=("$sam/$item")
                                done; cd ../../../
                        done

                        for elem in "${R2[@]}"; do str2="$str2,$elem"; done
                        str2=${str2#?}
		
		##Construct singles input string
                        str3=""
			R3=()
			for sam in $qryfolder/*; do cd $sam
                                for item in ./*"singles"."fasta"; do [ -f "$item" ] || continue
                                item=$(basename ${item}); R3+=("$sam/$item")
                                done; cd ../../../
                        done

                        for elem in "${R3[@]}"; do str3="$str3,$elem"; done
                        str3=${str3#?}

      ##Make a run bowtie script
      bowtie2 -f -x $dbfolder/$sub/$sub --no-unal -p $thread --local -1 "$str1" -2 "$str2" -U "$str3" -S $tgtfolder/$sub/bowtie.out.$sub.sam
      error_check "Bowtie 2 failed at lined 99!" 

	         elif [ "$mode" == "Single-end" ]; then
		##Construct R1 input string
                        str1=""
			R1=()
			for sam in $qryfolder/*; do cd $sam
                                for item in ./*"R1"."fasta"; do [ -f "$item" ] || continue
                                item=$(basename ${item}); R1+=("$sam/$item")
                                done; cd ../../../
                        done

                        for elem in "${R1[@]}"; do str1="$str1,$elem"; done
                        str1=${str1#?}
						
                        bowtie2 -f -x $dbfolder/$sub/$sub --un-unal -p $thread --local -U "$str1" -S $tgtfolder/$sub/bowtie.out.$sub.sam
                else
                        echo "No read mode selected!" >> 00.script/01.log/job.monitor_$sub.txt
                        exit
        fi						
#############################
	srcfolder="03.blast/03.bowtie.nucl/run.$run"
	tgtfolder="04.retrieve.reads/03.bowtie.nucl/run.$run"
	seqtype="nucl"
	logfolder="bowtie.log/bowtie.run.$run"
	sleeptime="20"
	unmap="map"

            mkdir -p $tgtfolder/$sub
			if [ "$mode" == "paired-end" ]; then
				perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.sam $unmap $tgtfolder/$sub/retrieved.$sub.R1.fasta $tgtfolder/$sub/unmap.$sub.R1.fasta $tgtfolder/$sub/retrieved.$sub.R2.fasta $tgtfolder/$sub/unmap.$sub.R2.fasta
				error_check "Bowtie broke at $sub"
			elif [ "$mode" == "single-end" ]; then                        
				perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.sam $unmap $tgtfolder/$sub/retrieved.$sub.fasta $tgtfolder/$sub/unmap.$sub.R1.fasta
			else
				echo "Error, specify read mode" >> 00.script/01.log/job.monitor_$sub.txt
				exit
            fi
fi

#Perform the Trinity
memory=24
tgtfolder=06.assembly/03.bowtie.nucl/run."$run"/"$sub".trinity
srcfolder=04.retrieve.reads/03.bowtie.nucl/run."$run"/"$sub"

while [[ ! -s $tgtfolder/Trinity.fasta ]]; do
		rm -r $tgtfolder		
		mkdir -p $tgtfolder

		Trinity  --seqType fa --max_memory "$memory"G --CPU 8 --left "$srcfolder"/retrieved."$sub".R1.fasta --right "$srcfolder"/retrieved."$sub".R2.fasta --output $tgtfolder --min_contig_length 100		
done

echo "Trinity output to $tgtfolder/$sub.trinity"

# Truncate the header
srcfolder="06.assembly/03.bowtie.nucl/run.$run"
time perl 00.script/06.truncate.header.pl $srcfolder/"$sub".trinity/Trinity.fasta $srcfolder/"$sub".trinity/Trinity.new.fasta

# Perform blastx
srcfolder="06.assembly/03.bowtie.nucl/run.$run"
dbfolder="01.data/05.SplitGenes/01.Protein/run.0"
tgtfolder="07.map.back/03.bowtie.nucl/run.$run"
reffolder="01.data/05.SplitGenes/01.Protein/run.0"
mkdir -p $tgtfolder/$sub
blastx -db $dbfolder/$sub/$sub.fasta -query $srcfolder/"$sub".trinity/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.out -evalue 1e-5  -outfmt 6 -num_threads 1 -max_target_seqs 1
blastx -db $dbfolder/$sub/$sub.fasta -query $srcfolder/"$sub".trinity/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.xml.out -evalue 1e-5  -outfmt 5 -num_threads 1 -max_target_seqs 1

# Perform blastn
srcfolder="06.assembly/03.bowtie.nucl/run.$run"
dbfolder="01.data/05.SplitGenes/02.Transcript/run.0"
tgtfolder="07.map.back/02.blastn/run.$run"
mkdir -p $tgtfolder/$sub
blastn -db $dbfolder/$sub/$sub.fasta -query $srcfolder/"$sub".trinity/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.out -evalue 1e-10  -outfmt 6 -num_threads 1 -max_target_seqs 1
blastn -db $dbfolder/$sub/$sub.fasta -query $srcfolder/"$sub".trinity/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.xml.out -evalue 1e-10  -outfmt 5 -num_threads 1 -max_target_seqs 1

## Script ends

#PBS -S /bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=48:HIGHMEM
#PBS -l walltime=00:45:00
#PBS -l mem=50gb
#PBS -N prerun
#PBS -o zzz.o%j
# this part only needs to be run once, it prepares and formats the data
# error-checking function, kills job and reports to error file on non-zero exit
error_check() {
	if [[ $? -ne 0 ]]; then
	echo "$1" >> job.monitor.txt
	exit $?
	fi
return 0
}

logfolder="00.script/01.log/"
#Load required modules
module load ncbiblast+
module load perl
module load mcl
module load R
module load diamond
module load python/2.7.8
# ncbi blast to create protein databases
makeblastdb -in 01.data/00.PriorData/proteome.fa -dbtype prot
error_check "Failed to make blast db, check preRun.sh line 31"

# protein blast
echo "Running blastp"
time blastp -num_threads 24 -db 01.data/00.PriorData/proteome.fa -query 01.data/00.PriorData/proteome.fa -out 01.data/03.MCL/01.blast/blast.all.out -evalue 1e-5 -outfmt 6
error_check "Failed protein blast, check preRun.sh line 36"

# mcl, cluster ortholog groups
# convert blast output into similarity data
echo "Preparing mcl graph"
time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.blast/blast.all.out 01.data/03.MCL/02.mcl/mcl.graph.txt
error_check "PERL failed us! Check preRun.sh line 42 and a1.mcl.prepare.graph.pl"

# mcl does the clustering to define gene family
echo "Performing gene family clustering"
time mcl 01.data/03.MCL/02.mcl/mcl.graph.txt --abc -o 01.data/03.MCL/02.mcl/mcl.out.txt -I 1.5
error_check "Gene family clustering failed, check preRun.sh line 48"

# construct meta-group, combine ortholog groups into meta-group, each group contains 1000 genes
echo "Running that one R script"
time Rscript 00.script/a3.geneSelection.R 01.data/03.MCL/02.mcl/mcl.out.txt 01.data/04.GeneOfInterest/GeneID.txt 1000 Potri
error_check "R failed us! Check preRun.sh line 54"

awk '{print} NR % 1000 == 0 {print ""}' 01.data/04.GeneOfInterest/GeneID.txt >> tmp && mv tmp 01.data/04.GeneOfInterest/GeneID.txt
error_check "AWK failed us! Check preRun.sh line 58"
# Split gene, based on the meta-group, Split gene sequences accordingly
echo "Running splitGene on protein"
time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/proteome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.SplitGenes/01.Protein/run.0 1000
error_check "PERL failed us! Check preRun.sh line 62 and a4.splitGene.pl"

echo "Running splitGene on transcript"
time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/transcriptome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.SplitGenes/02.Transcript/run.0 1000
error_check "PERL failed us! Check preRun.sh line 67 and a4.splitGene.pl"

# get meta-data for the meta-group, eg. gene/protein length, which group each gene belongs to
echo "getting relevant info"
time perl 00.script/a5.relevantInfo.pl 01.data/04.GeneOfInterest/GeneID.txt 01.data/00.PriorData/proteome.fa 01.data/00.PriorData/transcriptome.fa 01.data/04.GeneOfInterest/GeneID.v1.txt 1000
error_check "PERL failed us! Check preRun.sh line 73 and a5.releventInfo.pl"

#########################################################
#####derived from 01.folder.fastaCombinePairedEnd.pl##### !!!MOVE TO PRERUN!!!
#########################################################

echo 'Separating paired-end and single reads...' >> $logfolder/job.monitor_preRun.txt
##Get and loop through samples
for sam in 01.data/01.Fastq/*; do
  if [ -d $sam ]; then
	sam=$(basename ${sam})
        echo 'Running paired end script for $sam' >> $logfolder/job.monitor_preRun.txt
        time python2.7 00.script/01.fastaCombinePairedEnd.py 01.data/01.Fastq/$sam/$sam.R1.fastq 01.data/01.Fastq/$sam/$sam.R2.fastq " "
  fi
done

echo 'Finished separating reads!' >> $logfolder/job.monitor_preRun.txt

#########################################################
#########derived from 01.fastq2fasta.folder.pl########### !!!MOVE TO PRERUN!!!
#########################################################

echo 'Running 01.fastq2fasta.folder.pl ....' >> $logfolder/job.monitor_preRun.txt
mode="paired-end"
srcfolder="01.data/01.Fastq"
tgtfolder="01.data/02.Fasta"

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

echo 'Finished converting fastq to fasta!' >> $logfolder/job.monitor_preRun.txt

#########################################################
#########derived from 01.folder.IDConverter.pl########### !!!MOVE TO PRERUN!!!
#########################################################

echo 'Running 01.folder.IDConverter.pl ....' >> $logfolder/job.monitor_preRun.txt

## read in parameters required by the script
srcfolder="01.data/02.Fasta"
mode="paired-end"

for db in 01.data/02.Fasta/*; do
  if [ -d "$db" ]; then
        db=$(basename ${db})
        if [ $mode == "paired-end" ]; then
                time cat $srcfolder/$db/$db.R1.fasta_pairs_R1.fasta | awk '{if(NR%2==0){print $_} else{print ">"NR/2+0.5}}' > $srcfolder/$db/$db.R1.fasta_simple.fasta
                time cat $srcfolder/$db/$db.R2.fasta_pairs_R2.fasta | awk '{if(NR%2==0){print $_} else{print ">"NR/2+0.5}}' > $srcfolder/$db/$db.R2.fasta_simple.fasta
                time cat $srcfolder/$db/$db.R1.fasta_singles.fasta | awk '{if(NR%2==0){print $_} else{print ">"NR/2+0.5}}' > $srcfolder/$db/$db.singles.fasta_simple.fasta
    elif [ $mode == "single-end" ]; then
                cat $srcfolder/$db/$db.fasta | awk '{if(NR%2==0){print $_} else{print ">"NR/2+0.5}}' > $srcfolder/$db/$db.simple.fasta
        else
                echo 'Error: no read mode available!';
                exit 2
        fi
  fi
done

echo 'Finished converting folder ID!' >> $logfolder/job.monitor_prerun.txt

:>01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta
:>01.data/05.SplitGenes/03.Full.Length/full.length.contigs.prot.fasta

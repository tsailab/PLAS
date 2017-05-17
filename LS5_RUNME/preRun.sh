#!/bin/bash
#SBATCH -J PLAS_Prerun
#SBATCH -o PLAS_Prerun.o%j
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal
#SBATCH -t 3:00:00

#SBATCH --mail-user=sjq28742@uga.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# this part only needs to be run once, it prepares and formats the data
# error-checking function, kills job and reports to error file on non-zero exit
error_check() {
	if [[ $? -ne 0 ]]; then
	echo "$1" >> job.monitor.txt
	exit $?
	fi
return 0
}

#Load required modules
module load blast
module load perl

# ncbi blast to create protein databases
makeblastdb -in 01.data/00.PriorData/proteome.fa -dbtype prot
error_check "Failed to make blast db, check preRun.sh line 28"

# protein blast
echo "Running blastp"
time blastp -num_threads 24 -db 01.data/00.PriorData/proteome.fa -query 01.data/00.PriorData/proteome.fa -out 01.data/03.MCL/01.blast/blast.all.out -evalue 1e-5 -outfmt 6
error_check "Failed protein blast, check preRun.sh line 34"

# mcl, cluster ortholog groups
# convert blast output into similarity data
echo "Preparing mcl graph"
time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.blast/blast.all.out 01.data/03.MCL/02.mcl/mcl.graph.txt
error_check "PERL failed us! Check preRun.sh line 40 and a1.mcl.prepare.graph.pl"
wait

# mcl does the clustering to define gene family
echo "Performing gene family clustering"
time $WORK/bin/mcl 01.data/03.MCL/02.mcl/mcl.graph.txt --abc -o 01.data/03.MCL/02.mcl/mcl.out.txt -I 1.5
error_check "Gene family clustering failed, check preRun.sh line 46"
wait

# construct meta-group, combine ortholog groups into meta-group, each group contains 1000 genes
module load Rstats 
echo "Running that one R script"
time Rscript 00.script/a3.geneSelection.R 01.data/03.MCL/02.mcl/mcl.out.txt 01.data/04.GeneOfInterest/GeneID.txt 1000 Potri
error_check "R failed us! Check preRun.sh line 53"
wait

# Split gene, based on the meta-group, Split gene sequences accordingly
echo "Running splitGene on protein"
time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/proteome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.SplitGenes/01.Protein/run.0 1000
error_check "PERL failed us! Check preRun.sh line 59 and a4.splitGene.pl"
wait

echo "Running splitGene on transcript"
time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/transcriptome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.SplitGenes/02.Transcript/run.0 1000
error_check "PERL failed us! Check preRun.sh line 64 and a4.splitGene.pl"
wait

# get meta-data for the meta-group, eg. gene/protein length, which group each gene belongs to
echo "getting relevant info"
time perl 00.script/a5.releventInfo.pl 01.data/04.GeneOfInterest/GeneID.txt 01.data/00.PriorData/proteome.fa 01.data/00.PriorData/transcriptome.fa 01.data/04.GeneOfInterest/GeneID.v1.txt 1000
error_check "PERL failed us! Check preRun.sh line 70 and a5.releventInfo.pl"
wait

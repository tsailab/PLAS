#!/bin/bash

platform="Zcluster"
mode="paired-end"

####################### this starts with proteome #######################
## this part only needs to be run once, it prepares and formats the data
# WU-BLAST, construct similarity matrix for ortholog group clustering
## make blast database
/usr/local/wublast/latest/xdformat -p -o 01.data/03.MCL/01.blast/database 01.data/00.PriorData/proteome.fa
wait
## protein blast
# possibly figure out a way about parallelization
/usr/local/wublast/latest/blastp 01.data/03.MCL/01.blast/database 01.data/00.PriorData/proteome.fa -o 01.data/03.MCL/01.blast/wu.blast.all.out -e 1e-5 -mformat 2 -cpus 4 -wordmask seg
wait

# mcl, cluster ortholog groups
# convert blast output into similarity data
time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.blast/wu.blast.all.out 01.data/03.MCL/02.mcl/mcl.graph.txt wu
wait
# mcl does the clustering to define gene family
time /usr/local/mcl/latest/bin/mcl 01.data/03.MCL/02.mcl/mcl.graph.txt --abc -o 01.data/03.MCL/02.mcl/mcl.out.txt -I 1.5
wait

# construct meta-group, combine ortholog groups into meta-group, each group contains 1000 genes
module load R/3.1.2    # this is for Sapelo
#Using R/3.3.0, as 3.1.2 is not available on zcluster
time /usr/local/R/3.3.0/bin/Rscript 00.script/a3.geneSelection.R "01.data/03.MCL/02.mcl/mcl.out.txt" "01.data/04.GeneOfInterest/GeneID.txt" 1000 "Potri"
wait

# Split gene, based on the meta-group, Split gene sequences accordingly
time perl 00.script/a4.SplitGene.pl 01.data/00.PriorData/proteome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.SplitGenes/01.Protein/run.0 1000
wait
time perl 00.script/a4.SplitGene.pl 01.data/00.PriorData/transcriptome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.SplitGenes/02.Transcript/run.0 1000
wait

# get meta-data for the meta-group, eg. gene/protein length, which group each gene belongs to
time perl 00.script/a5.releventInfo.pl 01.data/04.GeneOfInterest/GeneID.txt 01.data/00.PriorData/proteome.fa 01.data/00.PriorData/transcriptome.fa 01.data/04.GeneOfInterest/GeneID.v1.txt 1000
wait

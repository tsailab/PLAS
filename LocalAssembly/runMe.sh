#!/bin/bash

platform="Zcluster"
mode="paired-end"

####################### this starts with proteome #######################
# WU-BLAST
#/usr/local/wublast/latest/xdformat -p -o 01.data/03.MCL/01.blast/database 01.data/00.PriorData/proteome.fa
#/usr/local/wublast/latest/blastp 01.data/03.MCL/01.blast/database 01.data/00.PriorData/proteome.fa -o 01.data/03.MCL/01.blast/blast.all.out -e 1e-5 -mformat 2 -cpus 4 -wordmask seg

# NCBI-BLAST
#/usr/local/ncbiblast+/latest/bin/makeblastdb -in 01.data/00.PriorData/proteome.fa prot
#time /usr/local/ncbiblast+/latest/bin/blastp -num_threads 4 -db 01.data/00.PriorData/proteome.fa -query 01.data/00.PriorData/proteome.fa -out 01.data/03.MCL/01.blast/blast.all.out -evalue 1e-5 -outfmt 6

# mcl
#time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.blast/ncbi.blast.all.out 01.data/03.MCL/02.mcl/mcl.graph.txt
#time /usr/local/mcl/latest/bin/mcl 01.data/03.MCL/02.mcl/mcl.graph.txt --abc -o 01.data/03.MCL/02.mcl/mcl.out.txt -I 1.5

#time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.blast/wu.blast.all.out 01.data/03.MCL/02.mcl/mcl.graph.txt wu
#time /usr/local/mcl/latest/bin/mcl 01.data/03.MCL/02.mcl/wu.mcl.graph.txt --abc -o 01.data/03.MCL/02.mcl/mcl.out.txt -I 1.5

# construct meta-group
#module load R/3.1.2
#time /usr/local/apps/R/3.1.2/bin/Rscript 00.script/a3.geneSelection.R 01.data/03.MCL/02.mcl/mcl.out.txt 01.data/04.GeneOfInterest/GeneID.txt 1000 Potri

# split gene
#time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/proteome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.splitGenes/01.Protein/run.0 1000
#time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/transcriptome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.splitGenes/02.Transcript/run.0 1000

#time perl 00.script/a5.releventInfo.pl 01.data/04.GeneOfInterest/GeneID.txt 01.data/00.PriorData/proteome.fa 01.data/00.PriorData/transcriptome.fa 01.data/00.PriorData/gene.gff3 01.data/04.GeneOfInterest/GeneID.v1.txt 1000

###### preprocess data

#time perl 00.script/02.makeblastdb.folder.pl 01.data/05.splitGenes/01.Protein/run.0 prot DIAMOND $platform
#time perl 00.script/02.makeblastdb.folder.pl 01.data/05.splitGenes/02.Transcript/run.0 nucl FALSE $platform
#time perl 00.script/01.folder.fastaCombinePairedEnd.pl  01.data/01.Fastq " " $platform
#time perl 00.script/01.fastq2fasta.folder.pl 01.data/01.Fastq 01.data/02.Fasta $mode $platform
#time perl 00.script/01.folder.IDConverter.pl 01.data/02.Fasta $mode $platform

mkdir -p 01.data/05.splitGenes/03.Full.Length/
:>01.data/05.splitGenes/03.Full.Length/full.length.contigs.nucl.fasta
:>01.data/05.splitGenes/03.Full.Length/full.length.contigs.prot.fasta

###################################################l
#: '
a=-1
evalue=1e-3
while [ $a -le 14 ]
do
    b=`expr $a + 1`
    echo "The run: $b" >> job.monitor.txt

    if [ $b -eq 0 ];then
        time perl 00.script/03.diamond.folder.pl 01.data/02.Fasta 01.data/05.splitGenes/01.Protein/run.$b 03.blast/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b $evalue $mode $platform 
        time perl 00.script/040.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.$b 04.retrieve.reads/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b 1000 genome $mode $platform 20
        echo "Have been here"
    else
        time perl 00.script/021.makebowtiedb.folder.pl 01.data/05.splitGenes/02.Transcript/run.$b $platform 10
        time perl 00.script/03.bowtie.folder.pl 01.data/02.Fasta 01.data/05.splitGenes/02.Transcript/run.$b 03.blast/03.bowtie.nucl/run.$b nucl local bowtie.log/bowtie.run.$b $mode $platform 10
        time perl 00.script/04.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.$b 04.retrieve.reads/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b $mode $platform 20
    fi

    time perl 00.script/06.assembly.trinity.folder.pl 04.retrieve.reads/03.bowtie.nucl/run.$b 06.assembly/03.bowtie.nucl/run.$b genome $mode $platform 20
    time perl 00.script/06.truncate.header.folder.pl 06.assembly/03.bowtie.nucl/run.$b $platform 20
    time perl 00.script/07.blastx.back.pl 06.assembly/03.bowtie.nucl/run.$b 01.data/05.splitGenes/01.Protein/run.0 07.map.back/03.bowtie.nucl/run.$b $platform 10
    time perl 00.script/07.blastn.back.pl 06.assembly/03.bowtie.nucl/run.$b 01.data/05.splitGenes/02.Transcript/run.0 07.map.back/02.blastn/run.$b $platform 10    

    c=`expr $b + 1`
    time perl 00.script/100.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.$b 07.map.back/03.bowtie.nucl/run.$b 01.data/05.splitGenes/02.Transcript/run.$c 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 
    time perl 00.script/10.folder.detect.full.length.seq.pl 07.map.back/03.bowtie.nucl/run.$b 06.assembly/03.bowtie.nucl/run.$b 01.data/05.splitGenes/01.Protein/run.0 01.data/05.splitGenes/03.Full.Length/run.$c pct 0.1
    cat 01.data/05.splitGenes/03.Full.Length/run.$c/full.length.contigs.nucl.fasta >> 01.data/05.splitGenes/03.Full.Length/full.length.contigs.nucl.fasta
    cat 01.data/05.splitGenes/03.Full.Length/run.$c/full.length.contigs.prot.fasta >> 01.data/05.splitGenes/03.Full.Length/full.length.contigs.prot.fasta
    a=`expr $a + 1`
done


#### summarize all runs
grep ">" 01.data/05.splitGenes/03.Full.Length/full.length.contigs.nucl.fasta > 01.data/05.splitGenes/03.Full.Length/count1
sed -i "s/>//" 01.data/05.splitGenes/03.Full.Length/count1
perl 00.script/b3.full.length.format.pl 01.data/04.GeneOfInterest/GeneID.v1.txt 01.data/05.splitGenes/03.Full.Length/count1 01.data/05.splitGenes/03.Full.Length/count2 01.data/05.splitGenes/03.Full.Length/count3
#'
#### assemble unmapped reads
mkdir -p 08.full.length
cp 01.data/05.splitGenes/03.Full.Length/count3 08.full.length/
cp 01.data/05.splitGenes/03.Full.Length/full.length.contigs.nucl.fasta 08.full.length/
time perl 00.script/c9.get.full.length.seq.pl 08.full.length/count3 08.full.length/full.length.contigs.nucl.fasta 08.full.length/Final.v1.fasta

: '
## start of second part
time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 32 -max_target_seqs 1

time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.v1.fasta -dbtype nucl
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.v1.fasta -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5

time perl 00.script/c11.remove.redundancy.pl 08.full.length/Final.v1.blastn.xml.out 08.full.length/Final.v1.fasta 08.full.length/Final.v2.fasta self 08.full.length/Final.v1.ptr.blastx.out > 08.full.length/record.run1

time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.v2.fasta -dbtype nucl
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.v2.fasta -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5
time perl 00.script/101.transfer.saturate.seq.pl 08.full.length/Final.v2.ptr.blastx.out 08.full.length/Final.v2.fasta 01.data/00.PriorData/ptr.proteome.fa 08.full.length Final.v2.ptr pct 0.02

cd 08.full.length/
ln -sf Final.v2.ptr.full.contigs.nucl.fasta Final.fasta
cd ../

#module load bowtie2/2.2.4
export PATH=$PATH:/usr/local/bowtie2/2.2.3/bin/
time bowtie2-build -f -q 08.full.length/Final.fasta 08.full.length/Final
time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 08.full.length/Final 09.bowtie.full.length unmap $platform
time perl 00.script/c10.unmapped.reads.trinity.pl 09.bowtie.full.length $platform
time perl 00.script/06.truncate.header.pl 10.unmapped.reads.trinity/Trinity.fasta 10.unmapped.reads.trinity/Trinity.new.fasta

## end of the second part

## make blast database for Mimulus and Ath
#time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 01.data/00.PriorData/ptr.proteome.fa -dbtype prot
#time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 01.data/00.PriorData/transcriptome.fa -dbtype nucl
#time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 01.data/00.PriorData/Ath.proteome.fa -dbtype prot
#time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 01.data/00.PriorData/uniprot_sprot.fasta -dbtype prot

### combine full length contigs from local assembly and assembly of unmapped reads
time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.fasta -dbtype nucl
#time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.fasta -query 10.unmapped.reads.trinity/Trinity.new.fasta -out 10.unmapped.reads.trinity/Trinity.new.blastn.out -evalue 1e-5 -outfmt 6 -num_threads 4 -max_target_seqs 1
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.fasta -query 10.unmapped.reads.trinity/Trinity.new.fasta -out 10.unmapped.reads.trinity/Trinity.new.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 1
time perl 00.script/c11.remove.redundancy.pl 10.unmapped.reads.trinity/Trinity.new.blastn.xml.out 10.unmapped.reads.trinity/Trinity.new.fasta 10.unmapped.reads.trinity/temp.fasta query > 10.unmapped.reads.trinity/remove.redundancy.log
mkdir 01.data/06.TargetTranscriptome
cat 08.full.length/Final.fasta 10.unmapped.reads.trinity/temp.fasta > 01.data/06.TargetTranscriptome/transcriptome.v1.fa
#rm 10.unmapped.reads.trinity/temp.fasta

time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 01.data/06.TargetTranscriptome/transcriptome.v1 14.evaluation/01.reads.mapping/01.Local map $platform
time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 05.Trinity/Trinity.new 14.evaluation/01.reads.mapping/02.Trinity map $platform

## evaluate assemblies
# blast to ptr proteome
time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 05.Trinity/Trinity.new.fasta -out 05.Trinity/Trinity.new.ptr.blastx.out -evalue 1e-5  -outfmt 6 -num_threads 32 -max_target_seqs 1
time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 01.data/06.TargetTranscriptome/transcriptome.v1.fa -out 01.data/06.TargetTranscriptome/transcriptome.v1.ptr.blastx.out -evalue 1e-5  -outfmt 6 -num_threads 32 -max_target_seqs 1

time perl 00.script/101.transfer.saturate.seq.pl 05.Trinity/Trinity.new.ptr.blastx.out 05.Trinity/Trinity.new.fasta 01.data/00.PriorData/ptr.proteome.fa 05.Trinity ptr pct 0.02
time perl 00.script/101.transfer.saturate.seq.pl 01.data/06.TargetTranscriptome/transcriptome.v1.ptr.blastx.out 01.data/06.TargetTranscriptome/transcriptome.v1.fa  01.data/00.PriorData/ptr.proteome.fa 01.data/06.TargetTranscriptome ptr pct 0.02

time perl 00.script/101.transfer.saturate.seq.pl 05.Trinity/Trinity.new.ath.blastx.out 05.Trinity/Trinity.new.fasta 01.data/00.PriorData/ath.proteome.fa 05.Trinity ath pct 0.02
time perl 00.script/101.transfer.saturate.seq.pl 01.data/06.TargetTranscriptome/transcriptome.v1.ath.blastx.out 01.data/06.TargetTranscriptome/transcriptome.v1.fa  01.data/00.PriorData/ath.proteome.fa 01.data/06.TargetTranscriptome ath pct 0.02

time perl 00.script/101.transfer.saturate.seq.pl 05.Trinity/Trinity.new.uniprot.blastx.out 05.Trinity/Trinity.new.fasta 01.data/00.PriorData/uniprot_sprot.fasta 05.Trinity uniprot pct 0.02
time perl 00.script/101.transfer.saturate.seq.pl 01.data/06.TargetTranscriptome/transcriptome.v1.uniprot.blastx.out 01.data/06.TargetTranscriptome/transcriptome.v1.fa  01.data/00.PriorData/uniprot_sprot.fasta 01.data/06.TargetTranscriptome uniprot pct 0.02

## compare local and trinity
# find unique recovered genes in both categories
time perl 00.script/c13.compare.Trinity.Local.pl 01.data/06.TargetTranscriptome/ptr.full.contigs.nucl.fasta 01.data/06.TargetTranscriptome/transcriptome.v1.ptr.blastx.out 05.Trinity/ptr.full.contigs.nucl.fasta 05.Trinity/Trinity.new.ptr.blastx.out 12.Local.Trinity.comparison/01.Protein 01.data/00.PriorData/ptr.proteome.fa
time perl 00.script/c13.compare.Trinity.Local.pl 01.data/06.TargetTranscriptome/ptr.full.contigs.nucl.fasta 01.data/06.TargetTranscriptome/transcriptome.v1.ptr.blastn.out 05.Trinity/ptr.full.contigs.nucl.fasta 05.Trinity/Trinity.new.ptr.blastn.out 12.Local.Trinity.comparison/02.Transcript 01.data/00.PriorData/ptr.transcriptome.fa

time perl 00.script/c14.retrieveseq.Trinity.Local.pl 12.Local.Trinity.comparison/Trinity_Local_results.txt 13.Local.Trinity.alignment 01.data/06.TargetTranscriptome/transcriptome.v1.fa 01.data/06.TargetTranscriptome/proteome.v1.fa 05.Trinity/Trinity.new.fasta 05.Trinity/Trinity.prot.fasta 01.data/00.PriorData/ptr.proteome.fa
time perl 00.script/c15.clustalo.Trinity.Local.pl 13.Local.Trinity.alignment/01.data 13.Local.Trinity.alignment/02.alignment

time perl 00.script/c16.find.hybridicity.pl 12.Local.Trinity.comparison/01.Protein/Trinity_Local_intersect.txt 12.Local.Trinity.comparison/01.Protein/Trinity_Local_intersect_quality.v2.txt
time perl 00.script/c16.find.hybridicity.pl 12.Local.Trinity.comparison/02.Transcript/Trinity_Local_intersect.txt 12.Local.Trinity.comparison/02.Transcript/Trinity_Local_intersect_quality.v2.txt

## evaluate assembly
# map reads
export PATH=$PATH:/usr/local/bowtie2/2.2.3/bin/
time bowtie2-build -f -q 01.data/06.TargetTranscriptome/transcriptome.v1.fa 01.data/06.TargetTranscriptome/transcriptome.v1
time bowtie2-build -f -q 05.Trinity/Trinity.new.fasta 05.Trinity/Trinity.new

time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 01.data/06.TargetTranscriptome/transcriptome.v1 14.evaluation/01.reads.mapping/01.Local map $platform
time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 05.Trinity/Trinity.new 14.evaluation/01.reads.mapping/02.Trinity map $platform

# comput Nx statistcs
mkdir -p 14.evaluation/02.NxStatistics
/usr/local/trinity/r20140717/util/TrinityStats.pl 01.data/06.TargetTranscriptome/transcriptome.v1.fa > 14.evaluation/02.NxStatistics/01.Local.Nx.log
/usr/local/trinity/r20140717/util/TrinityStats.pl  05.Trinity/Trinity.new.fasta > 14.evaluation/02.NxStatistics/02.Trinity.log

/usr/local/trinity/r20140717/util/TrinityStats.pl 01.data/06.TargetTranscriptome/ptr.full.contigs.nucl.fasta > 14.evaluation/02.NxStatistics/01.Local.full.log
/usr/local/trinity/r20140717/util/TrinityStats.pl 05.Trinity/ptr.full.contigs.nucl.fasta > 14.evaluation/02.NxStatistics/02.Trinity.full.log

/usr/local/trinity/r20140717/util/TrinityStats.pl 01.data/06.TargetTranscriptome/ptr.incomplete.contigs.nucl.fasta > 14.evaluation/02.NxStatistics/01.Local.incomplete.log
/usr/local/trinity/r20140717/util/TrinityStats.pl 05.Trinity/ptr.incomplete.contigs.nucl.fasta > 14.evaluation/02.NxStatistics/02.Trinity.incomplete.log

/usr/local/trinity/r20140717/util/TrinityStats.pl 01.data/06.TargetTranscriptome/ptr.unmapped.contigs.nucl.fasta > 14.evaluation/02.NxStatistics/01.Local.unmapped.log
/usr/local/trinity/r20140717/util/TrinityStats.pl 05.Trinity/ptr.unmapped.contigs.nucl.fasta > 14.evaluation/02.NxStatistics/02.Trinity.unmapped.log

# estimate transcript abundance
export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
export PATH=/usr/local/gmap-gsnap/latest/bin/:${PATH}
export PATH=${PATH}:/usr/local/cd-hit/4.6.1-2012-08-27/:/usr/local/rsem/1.2.20/:/usr/local/express/1.5.1/
time /usr/local/trinity/r20140717/util/align_and_estimate_abundance.pl --transcripts 01.data/06.TargetTranscriptome/transcriptome.v1.fa --est_method eXpress --aln_method bowtie2 --prep_reference
time /usr/local/trinity/r20140717/util/align_and_estimate_abundance.pl --transcripts 05.Trinity/Trinity.new.fasta --est_method eXpress --aln_method bowtie2 --prep_reference
time perl 00.script/13.estimate.abundance.pl 01.data/01.Fastq 14.evaluation/03.abundance/01.Local 01.data/06.TargetTranscriptome/transcriptome.v1.fa eXpress paired-end zcluster
time perl 00.script/13.estimate.abundance.pl 01.data/01.Fastq 14.evaluation/03.abundance/02.Trinity 05.Trinity/Trinity.new.fasta eXpress paired-end zcluster
/usr/local/trinity/2.0.6/util/abundance_estimates_to_matrix.pl --est_method eXpress 14.evaluation/03.abundance/01.Local/1T/results.xprs 14.evaluation/03.abundance/01.Local/2T/results.xprs 14.evaluation/03.abundance/01.Local/3T/results.xprs --name_sample_by_basedir --out_prefix 14.evaluation/03.abundance/01.Local/matrix
/usr/local/trinity/2.0.6/util/abundance_estimates_to_matrix.pl --est_method eXpress 14.evaluation/03.abundance/02.Trinity/1T/results.xprs 14.evaluation/03.abundance/02.Trinity/2T/results.xprs 14.evaluation/03.abundance/02.Trinity/3T/results.xprs --name_sample_by_basedir --out_prefix 14.evaluation/03.abundance/02.Trinity/matrix

perl 00.script/contig_ExN50_statistic.pl 14.evaluation/03.abundance/01.Local/matrix.TMM.fpkm.matrix 01.data/06.TargetTranscriptome/transcriptome.v1.fa | tee ExN50.stats > 14.evaluation/02.NxStatistics/03.Local.ExN50.log
perl 00.script/contig_ExN50_statistic.pl 14.evaluation/03.abundance/02.Trinity/matrix.TMM.fpkm.matrix 05.Trinity/Trinity.new.fasta | tee ExN50.stats > 14.evaluation/02.NxStatistics/04.Trinity.ExN50.log
'

## annotation
mkdir -p  16.annotation
#time perl 00.script/d18.annotation.pl 01.data/06.TargetTranscriptome/transcriptome.v1.fa 16.annotation/annotation.txt 01.data/06.TargetTranscriptome/transcriptome.v1.ptr.blastn.out 01.data/00.PriorData/Ptrichocarpa_210_v3.0.annotation_info.txt 01.data/06.TargetTranscriptome/transcriptome.v1.ath.blastx.out 01.data/00.PriorData/Athaliana_167_TAIR10.annotation_info.txt 01.data/06.TargetTranscriptome/transcriptome.v1.uniprot.blastx.out 01.data/00.PriorData/uniprot_sprot.fasta
time perl 00.script/d18.annotation.pl 01.data/06.TargetTranscriptome/transcriptome.v1.fa 16.annotation/annotation.txt 01.data/06.TargetTranscriptome/transcriptome.v1.ptr.blastx.out 01.data/00.PriorData/Ptrichocarpa_210_v3.0.annotation_info.txt 01.data/06.TargetTranscriptome/transcriptome.v1.ath.blastx.out 01.data/00.PriorData/Athaliana_167_TAIR10.annotation_info.txt 01.data/06.TargetTranscriptome/transcriptome.v1.uniprot.blastx.out 01.data/00.PriorData/uniprot_sprot.fasta

### estimate transcript abundance
export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
export PATH=/usr/local/gmap-gsnap/latest/bin/:${PATH}
time /usr/local/trinity/r20140717/util/align_and_estimate_abundance.pl --transcripts 01.data/06.TargetTranscriptome/transcriptome.v1.fa --est_method RSEM --aln_method bowtie2 --prep_reference
time perl 00.script/13.estimate.abundance.pl 01.data/01.Fastq 13.abundance/genome/RSEM 01.data/06.TargetTranscriptome/transcriptome.v1.fa RSEM paired-end zcluster
time perl 00.script/13.estimate.abundance.pl 01.data/01.Fastq 13.abundance/genome/eXpress 01.data/06.TargetTranscriptome/transcriptome.v1.fa eXpress paired-end zcluster
/usr/local/R/3.2.3/bin/Rscript 00.script/a8.summarize.expression.R 13.abundance/genome/RSEM 13.abundance/genome/RSEM/gene.fpkm.rsem.txt 01.data/00.PriorData/tissue_record.txt RSEM
/usr/local/R/3.2.3/bin/Rscript 00.script/a8.summarize.expression.R 13.abundance/genome/eXpress 13.abundance/genome/eXpress/gene.fpkm.express.txt 01.data/00.PriorData/tissue_record.txt eXpress



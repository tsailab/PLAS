#!/bin/bash

platform="Zcluster"
mode="paired-end"

###### preprocess data

## manually run
# construct DIAMOND and blast database for the meta-group
time perl 00.script/02.makeblastdb.folder.pl 01.data/05.splitGenes/01.Protein/run.0 prot DIAMOND $platform
wait
time perl 00.script/02.makeblastdb.folder.pl 01.data/05.splitGenes/02.Transcript/run.0 nucl FALSE $platform
wait
# separated paired reads and single reads
time perl 00.script/01.folder.fastaCombinePairedEnd.pl  01.data/01.Fastq " " $platform
wait
# convert fastq file to fasta
time perl 00.script/01.fastq2fasta.folder.pl 01.data/01.Fastq 01.data/02.Fasta $mode $platform
wait
# simplify read IDs
time perl 00.script/01.folder.IDConverter.pl 01.data/02.Fasta $mode $platform
wait
## automatic run
mkdir -p 01.data/05.splitGenes/03.Full.Length/
:>01.data/05.splitGenes/03.Full.Length/full.length.contigs.nucl.fasta
:>01.data/05.splitGenes/03.Full.Length/full.length.contigs.prot.fasta

###################################################l
#: '
a=0
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

## everything below needs to be run manually
time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 32 -max_target_seqs 1

time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.v1.fasta -dbtype nucl
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.v1.fasta -query 08.full.length/Final.v1.fasta -out 08.full.length/Final.v1.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5
time perl 00.script/c11.remove.redundancy.pl 08.full.length/Final.v1.blastn.xml.out 08.full.length/Final.v1.fasta 08.full.length/Final.v2.fasta self 08.full.length/Final.v1.ptr.blastx.out > 08.full.length/record.run1

time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.v2.fasta -dbtype nucl
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.v2.fasta -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 5
time /usr/local/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ptr.proteome.fa -query 08.full.length/Final.v2.fasta -out 08.full.length/Final.v2.ptr.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 32 -max_target_seqs 1
time perl 00.script/101.transfer.saturate.seq.pl 08.full.length/Final.v2.ptr.blastx.out 08.full.length/Final.v2.fasta 01.data/00.PriorData/ptr.proteome.fa 08.full.length Final.v2.ptr pct 0.02
cd 08.full.length/
ln -sf Final.v2.ptr.full.contigs.nucl.fasta Final.fasta
cd ../

#module load bowtie2/2.2.4
export PATH=$PATH:/usr/local/bowtie2/2.2.3/bin/
time bowtie2-build -f -q 08.full.length/Final.fasta 08.full.length/Final
wait
time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 08.full.length/Final 09.bowtie.full.length unmap $platform
wait
time perl 00.script/c10.unmapped.reads.trinity.pl 09.bowtie.full.length $platform
wait
time perl 00.script/06.truncate.header.pl 10.unmapped.reads.trinity/Trinity.fasta 10.unmapped.reads.trinity/Trinity.new.fasta


### combine full length contigs from local assembly and assembly of unmapped reads
time /usr/local/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.fasta -dbtype nucl
time /usr/local/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.fasta -query 10.unmapped.reads.trinity/Trinity.new.fasta -out 10.unmapped.reads.trinity/Trinity.new.blastn.xml.out -evalue 1e-5 -outfmt 5 -max_target_seqs 1
time perl 00.script/c11.remove.redundancy.pl 10.unmapped.reads.trinity/Trinity.new.blastn.xml.out 10.unmapped.reads.trinity/Trinity.new.fasta 10.unmapped.reads.trinity/temp.fasta query > 10.unmapped.reads.trinity/remove.redundancy.log
mkdir 01.data/06.TargetTranscriptome
cat 08.full.length/Final.fasta 10.unmapped.reads.trinity/temp.fasta > 01.data/06.TargetTranscriptome/transcriptome.v1.fa
#rm 10.unmapped.reads.trinity/temp.fasta


#PBS -S /bin/bash
#PBS -N assembleMe
#PBS -q highmem_q
#PBS -l nodes=1:ppn=8
#PBS -l mem=24gb
#PBS -l walltime=04:00:00
#PBS -d ./

run=${RUN}
let c=$run+1
cd $PBS_O_WORKDIR


#Load required modules
module load Perl/5.20.3-foss-2016b
module load BioPerl/1.7.1-foss-2016b-Perl-5.24.1
mkdir 01.data/05.SplitGenes/03.Full.Length/run.$c

#Transfer saturate seq
time perl 00.script/100.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.$run 07.map.back/03.bowtie.nucl/run.$run 01.data/05.SplitGenes/02.Transcript/run.$c 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 

#Detect full-length transcripts
time perl 00.script/10.folder.detect.full.length.seq.pl 07.map.back/03.bowtie.nucl/run.$run 06.assembly/03.bowtie.nucl/run.$run 01.data/05.SplitGenes/01.Protein/run.0 01.data/05.SplitGenes/03.Full.Length/run.$c pct 0.1

cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.nucl.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.nucl.fasta
cat 01.data/05.SplitGenes/03.Full.Length/run.$c/full.length.contigs.prot.fasta >> 01.data/05.SplitGenes/03.Full.Length/full.length.contigs.prot.fasta

## Script ends
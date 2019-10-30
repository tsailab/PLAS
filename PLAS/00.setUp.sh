#Make nesensary folders for PLAS 
mkdir 00.script/
mkdir 00.script/01.log
mkdir -p 01.data/{00.PriorData,01.Fastq,03.MCL,04.GeneOfInterest,05.SplitGenes,06.TargetTranscriptome}
mkdir 01.data/03.MCL/{01.blast,02.mcl}
mkdir 01.data/05.SplitGenes/{01.Protein,02.Transcript,03.Full.Length}
mkdir 01.data/05.SplitGenes/01.Protein/run.0
mkdir 01.data/05.SplitGenes/02.Transcript/run.0
mkdir -p 03.blast/03.bowtie.nucl 04.retrieve.reads/03.bowtie.nucl 06.assembly/03.bowtie.nucl 07.map.back/{02.blastn,03.bowtie.nucl}
mkdir -p 08.full.length 09.fulllength.bowtie 10.unmappedreads.trinity



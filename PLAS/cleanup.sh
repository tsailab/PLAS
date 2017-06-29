#!/bin/bash

mv job.monitor.txt archive/
mv PLAS_Prerun*.o* archive/
mv PLAS_runMe*.o* archive/
rm checkit.txt
echo "Removing leftover files from 03.blast/03.bowtie.nucl..."
rm -r 03.blast/03.bowtie.nucl/run.1
rm -r 03.blast/03.bowtie.nucl/run.2
rm -r 03.blast/03.bowtie.nucl/run.3

echo "Removing leftover files from 01.data/05.SplitGenes/02.Transcript..."
rm -r 01.data/05.SplitGenes/02.Transcript/run.1
rm -r 01.data/05.SplitGenes/02.Transcript/run.2
rm -r 01.data/05.SplitGenes/02.Transcript/run.3
rm -r 01.data/05.SplitGenes/02.Transcript/run.4

echo "Removing leftover files from 01.data/05.SplitGenes/03.Full.Length..."
rm -r 01.data/05.SplitGenes/03.Full.Length/run.1
rm -r 01.data/05.SplitGenes/03.Full.Length/run.2
rm -r 01.data/05.SplitGenes/03.Full.Length/run.3
rm -r 01.data/05.SplitGenes/03.Full.Length/run.4

echo "Removing leftover files from 04.retrieve.reads/03.bowtie.nucl..."
rm -r 04.retrieve.reads/03.bowtie.nucl/run.1
rm -r 04.retrieve.reads/03.bowtie.nucl/run.2
rm -r 04.retrieve.reads/03.bowtie.nucl/run.3

echo "Removing leftover files from 06.assembly/03.bowtie.nucl..."
rm -r 06.assembly/03.bowtie.nucl/run.1
rm -r 06.assembly/03.bowtie.nucl/run.2
rm -r 06.assembly/03.bowtie.nucl/run.3

echo "Removing leftover files from 00.script/10.transfer.script..."
rm -r 00.script/10.transfer.script/run.1
rm -r 00.script/10.transfer.script/run.2
rm -r 00.script/10.transfer.script/run.3

echo "Removing leftover files from 07.map.back/02.blastn..."
rm -r 07.map.back/02.blastn/run.1
rm -r 07.map.back/02.blastn/run.2
rm -r 07.map.back/02.blastn/run.3

echo "Removing leftover files from 07.map.back/03.bowtie.nucl..."
rm -r 07.map.back/03.bowtie.nucl/run.1
rm -r 07.map.back/03.bowtie.nucl/run.2
rm -r 07.map.back/03.bowtie.nucl/run.3

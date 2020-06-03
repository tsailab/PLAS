# Programmatic Structure of the PLAS Pipeline

## Overview

This series of scripts involves subscripts (often Perl), intermediate files, and scripts that call other scripts.

As such, this document is intended to clarify each programmatic step, parameter, dependency, and output of each script
and subscript, beginning with the setup script.  

There is an "all in one" bash script that involves blocking waits with `afterok` to run the scripts in order 
`01.PLAS.sh`.

Script files are located in the plain `PLAS` directory and the Perl, Python, and R subscripts called by the main bash
scripts are located in the `00.script` directory.

## 00 setUp.sh

**Summary**

Creates the directories and subdirectories used in the pipeline.

**Creates**

Creates following directory structure (there are some inconsistent numbers and it also looks like the used the same 
folder mames if they were shared by a process, likely due to a calling Perl script):

```bash
00.script
01.data
|_00.PriorData
|_01.Fastq
|_03.MCL
    |_01.blast
    |_02.mcl
|_04.GeneOfInterest
|_05.SplitGenes
    |_01.Protein
        |_run.0
    |_02.Transcript
        |_run.0
    |_03.Full.Length
|_06.TargetTranscriptome
03.blast
|_03.bowtie.nucl
04.retrieve.reads
|_03.bowtie.nucl
06.assembly
|_03.bowtie.nucl
07.map.back
|_02.blastn
|_03.bowtie.nucl
08.full.length
09.fulllength.bowtie
10.unmappedreads.trinity
```

**Calls**

N/A

## 01 PLAS.sh

This is the master script that calls the other scripts using the Sapelo2 blocking wait `afterok`, which waits for a
process to terminate before proceeding.

## 02.preRun.sh


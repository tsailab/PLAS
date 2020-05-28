# PLAS (ch 2) - Parallelized Local De Novo Assembly of Sequences

- [PLAS (ch 2) - Parallelized Local De Novo Assembly of Sequences](#plas-ch-2---parallelized-local-de-novo-assembly-of-sequences)
	- [Overview](#overview)
	- [Methods](#methods)
		- [Overview](#overview-1)
		- [Reference Organization](#reference-organization)
		- [Assembly of Conserved Sequences](#assembly-of-conserved-sequences)
		- [Assembly of Species Specific and Partial Sequences](#assembly-of-species-specific-and-partial-sequences)
		- [Assembly Evaluation](#assembly-evaluation)
	- [Datasets](#datasets)

## Overview

Reference based mapping and de novo (starting anew) assembly imrpoves computing efficiency and assembly quality

Input: filtered RNA-seq reads anda  reference proteome from a closely-related species as input

Reference proteome is first clustered by gene family before read mapping to group the input RNA seq data into bins for local de novo assembly

## Methods

### Overview

Input:

1. RNA seq read data in fastq format
2. Reference proteome/transcriptome from a closely-related species in fasta format

2 major pipeline components:

1. Assembly of conserved (mapped) sequences to full length
2. Assembles diverged (unmapped) sequences, species specific sequences, and partial sequences

Final Output: 

* De Novo assembled transcriptome

### Reference Organization

Markov Cluster Algorithm (`MCL Classification`) defines the organization of ther reference transcriptome into orthologous gene families

* Similarity measurement is the -log10 E value from the results of `WU-BLAST 2.2.6`

Highly similar genes must be be classified into the same group to avoid cross-assembly

OrthoMCL sorted gene families are combined into further meta-groups that results in a roughly equal number of genes &#8594; de novo assembly is performed on the read set by bin aligned to the meta group

### Assembly of Conserved Sequences

`DIAMOND` maps the input reads to there reference proteome.  They are then organized into bins based on the prescence/absence of significant huts against each meta group.  E value 1e-3 defines the sigificant cutoff value.

`Trinity version r020140717` assembles the independent read bins in parallel

Resulting contigs are mapped to the reference sequence with `BLAST` and only contigs with significant hits (E value 1e-5) are retained

Assembled contigs are the starting point for the next round of assembly.  `Bowtie2` performs read mapping and rebinned reads are used for the second iteration of de novo assembly

Process repeats depending on user specifications.

Results contigs (in parallel bins) are `BLASTN-mapped` against each other to remove redundency

* Sharing an overall identity >= 95% or when 90% of shorter sequences align with >= 90% identity across alignment = redundant

Product: intermediate assembly II - conserved stranscripts recovered to full length

* Full length transcript = aligned portion is >= 98% of the reference coding sequence of the same species or >= 90% of the protein reference of a different species

### Assembly of Species Specific and Partial Sequences

`Bowtie2` assembles divergent, species specific, and partial sequences.  Input reads are mapped to Intermediate Assembly II and only unmapped reads are retained

`Trinity` assembles unmapped reads de novo for Intermediate Assembly II

Intermediate Assembly III is then `BLASTN-searched` against the Intermedate Asssembly II to remove redundency

The combined non-redundant set is the final PLAS assembled transcriptome

### Assembly Evaluation

1. Assembly score and associated components (defined by `TransRate v1.0.1`)
2. The number of fully recovered transcripts and corresponding gene models
3. The alignment quality between the assemblies and the reference
4. The reconstruction of highly similar genes


## Datasets

Populus trichocarpa

* 3 replicates 21M paired end 75bp (PE75)

Cornus florida 'Appalachian Spring' (dogwood)

* 4 replicates 72M PE75

Lagerstroemia sp. (crepe myrtle)

* 3 replicates 94M PE75

Quercus robur (oak)

* 2 repliceates 19M PE75

Prunus persica (peach)

* 2 repliceates 33M PE75 

Prunus mume (plum) 

* 2 replicates 13M PE75

Salix purpurea L. (willow)

* 2 replicates 35M PE75


Populus trichocarpa v Arabidopsis


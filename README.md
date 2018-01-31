# Project directory
PLAS: Parallelized Local de novo Assembly of Sequences
PLAS is a pipeline for reference-guided de novo assembly of RNAseq data that boasts improved performance compared to the Trinity and CLC assembly pipelines.
PLAS leverages the conservation inherent to protein sequences as a reference, allowing the organization of input reads into independent bins before assembly is performed.
This pre-organization of reads simplifies the assembly process by reducing data complexity and facilitating parallelization, allowing PLAS to achieve higher coverage, accuracy, and computational efficiency.

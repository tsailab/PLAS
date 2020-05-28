# Programmatic Structure of the PLAS Pipeline

## Overview

This series of scripts involves subscripts (often Perl), intermediate files, and scripts that call other scripts.

As such, this document is intended to clarify each programmatic step, parameter, dependency, and output of each script
and subscript, beginning with the setup script.  

There is an "all in one" bash script that involves blocking waits to run the scripts in order, however, for the purpose
of this document the scripts will be treated as though they must be manually run from the outer level.  At no time
should a Perl subscript ever be manually called despite their input parameters listed.  These subscripts perform
complicated and often undocumented bioinformatics tasks and their steps will attempt be documented.

Script files are located in the plain `PLAS` directory and the Perl, Python, and R subscripts called by the main bash
scripts are located in the `00.script` directory.


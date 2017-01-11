#!/bin/bash

folder="06.assembly/03.bowtie.nucl"
for r in $folder/*
do
    for f in $r/*
    do
        echo $f
        rm $f/*.fa
        rm $f/*.finished
        rm $f/*_count
        rm $f/*bam*
        rm $f/*.txt
        rm $f/*.fa.*
        rm $f/*sam
        rm $f/*ebwt
        rm $f/*thread*
        rm -r $f/chrysalis
    done
done

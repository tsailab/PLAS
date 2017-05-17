#!/bin/bash
# Use -gt 1 to consume two arguments per pass in the loop (e.g. each
# argument has a corresponding value to go with it).
# Use -gt 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).
# note: if this is set to -gt 0 the /etc/hosts part is not recognized ( may be a bug )
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -d|--datadir)
    DATADIR="$2"
    shift # past argument
    ;;
    -w|--workingdir)
    WORKINGDIR="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done # while loop

data=${DATADIR}
wd=${WORKINGDIR}
cd $wd
master="Shell_master_QC.sh"
printf "#"'!'/bin/bash"\n" >$master
index=0
for file in $data/*.fq
do
  shortname=${file##*/}
  sample=${shortname%%.fq}
  out=$sample
  index=$(($index+1))
  sh_worker="run"$index"_"$out"_run.sh"
  printf "qsub "$sh_worker"\n" >>$master
  printf "#"'!'/bin/bash"\n" >$sh_worker
  printf "#PBS -N fqQC\n" >>$sh_worker
  printf "#PBS -q batch\n" >>$sh_worker 
  printf "#PBS -l nodes=1:ppn=8:HIGHMEM\n" >>$sh_worker 
  printf "#PBS -l walltime=10:00:00\n" >>$sh_worker 
  printf "cd "$wd"\n">>$sh_worker 
  printf "module load java/jdk1.8.0_20 fastqc\n" >>$sh_worker  
  printf "time fastqc  "$file" -t 8 -o   ./ ">>$sh_worker  
done


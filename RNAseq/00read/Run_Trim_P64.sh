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




### after getting the input


data=${DATADIR}
wd=${WORKINGDIR}
cd $wd
master="Shell_master_Trim64.sh"
printf "#"'!'/bin/bash"\n" >$master
index=0
for file in $data/*_1.fq
do
  sample=${file%%_1.fq}
  shortname=${sample##*/}
  out=$shortname
  index=$(($index+1))
  read1=$file
  read2=$sample"_2.fq" 
  r1_pair=$out"_P_1.fq.gz"
  r1_unpair=$out"_S_1.fq.gz"
  r2_pair=$out"_P_2.fq.gz"
  r2_unpair=$out"_S_2.fq.gz"  
  sh_worker="run"$index"_"$out"_run.sh"
  printf "qsub "$sh_worker"\n" >>$master
  printf "#"'!'/bin/bash"\n" >$sh_worker
  printf "#PBS -N fqTrim"$index"\n" >>$sh_worker
  printf "#PBS -q batch\n" >>$sh_worker 
  printf "#PBS -l nodes=1:ppn=10:HIGHMEM\n" >>$sh_worker 
  printf "#PBS -l walltime=20:00:00\n" >>$sh_worker 
  printf "cd "$wd"\n">>$sh_worker 
  printf "module load java/jdk1.8.0_20\n" >>$sh_worker  
  printf "time java -jar /usr/local/apps/trimmomatic/0.33/trimmomatic-0.33.jar  PE -threads 10 -phred64 \\">>$sh_worker  
  printf "\n">>$sh_worker  
  printf $read1" "$read2" "$r1_pair" "$r1_unpair" "$r2_pair" "$r2_unpair" \\">>$sh_worker  
  printf "\n">>$sh_worker  
  printf "ILLUMINACLIP:/usr/local/apps/trimmomatic/latest/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 \n">>$sh_worker  
  printf "zcat "$r1_pair" "$r1_unpair" >"$out"_1_CA.fastq\n">>$sh_worker 
  printf "zcat "$r2_pair" "$r2_unpair" >"$out"_2_CA.fastq\n">>$sh_worker 
done


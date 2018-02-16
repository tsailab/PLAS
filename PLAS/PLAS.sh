#PBS -S /bin/bash
#PBS -q batch
#PBS -N Master_PLAS
#PBS -l nodes=1:ppn=04:HIGHMEM
#PBS -l walltime=08:00:00
#PBS -l mem=50gb
cd $PBS_O_WORKDIR
repeats=15
counter=0
WAIT1=""
WAIT2=""
WAIT3=""
WAIT4=""

./preRun.sh

for sub in 01.data/05.SplitGenes/01.Protein/run.0/*; do 
	sub=$(basename ${sub})
	WAIT1+=$(qsub -N PLAS_$sub -v SUB=$sub, runMe.sh)
	WAIT1+=",afterok:"
done

WAIT1=${WAIT1%?????????}

while [ $repeats -ge $counter ]; do
	WAIT2=""
	for sub in 01.data/05.SplitGenes/01.Protein/run.0/*; do
		sub=$(basename ${sub})
		WAIT2+=$(qsub -N PLASL_${counter}_$sub -W depend=afterok:$WAIT1 -v RUN=$counter,SUB=$sub loopMe.sh)
		WAIT2+=",afterok:"
	done

WAIT2=${WAIT2%?????????}	

WAIT1=$(qsub -W depend=afterok:$WAIT2 -N PLASA_${counter} -v RUN=$counter, assembleMe.sh)
	let counter=$counter+1
done

WAIT3=$(qsub -W depend=afterok:$WAIT1 -N PLAS_mapped mapped.sh)

for dir in 01.data/01.Fastq/*; do
	if [ -d $dir ]; then
	dir=$(basename $dir)
		WAIT4+=$(qsub -W depend=afterok:$WAIT3 -N PLAS_unmapped -v SUB=$dir unmappedBowtie.sh)	
		WAIT4+=",afterok:"
	fi
done
WAIT4=${WAIT4%?????????}

qsub -W depend=afterok:$WAIT4 -N PLAS_final final.sh
	
echo "All done!"

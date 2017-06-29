#!/bin/bash
#SBATCH -J PLAS_Master_run
#SBATCH -o PLAS_Master_run.o%j
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal
#SBATCH -t 12:00:00

#SBATCH --mail-user=sjq28742@uga.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
echo "Running folder setUp" >> job.monitor.txt
./setUp.sh
echo "Running preRun" >> job.monitor.txt
./preRun.sh
echo "Running stage 1" >> job.monitor.txt
./1.runMe.sh
echo "Running stage 2" >> job.monitor.txt
./2.runMe.sh
echo "Running stage 3" >> job.monitor.txt
./3.runMe.sh
echo "Running stage 4" >> job.monitor.txt
./4.runMe.sh
echo "Running stage 5" >> job.monitor.txt
./5.runMe.sh
echo "Running stage 6" >> job.monitor.txt
./6.runMe.sh

#!/bin/bash
#SBATCH --job-name=geomx
#SBATCH --nodes=1
#SBATCH --nodelist=gpu03
##SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G
#SBATCH --partition=gpgpuq
#SBATCH --output=output

pwd; hostname; date
SECONDS=0

Rscript $1;
#python $1;

duration=$SECONDS
printf "\n$((duration / 60))m $((duration % 60))s ELAPSED\n"

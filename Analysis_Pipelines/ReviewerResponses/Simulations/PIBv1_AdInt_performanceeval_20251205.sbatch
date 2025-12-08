#!/bin/bash
#SBATCH -J AdIntSimEval
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96000
#SBATCH -t 24:00:00
#SBATCH -e PIBv1_AdInt_performanceeval_20251205_%j.stderr
#SBATCH -o PIBv1_AdInt_performanceeval_20251205_%j.stdout
#SBATCH --mail-user patrick.f.reilly@yale.edu
#SBATCH --mail-type ALL
#SBATCH -p day

#For mccleary:
module load R/4.3.2-foss-2022b

pushd /gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes/PIBv1_ReviewerResponses/Simulations/

/usr/bin/time -v ./PIBv1_AdInt_performanceeval_20251205.R

popd

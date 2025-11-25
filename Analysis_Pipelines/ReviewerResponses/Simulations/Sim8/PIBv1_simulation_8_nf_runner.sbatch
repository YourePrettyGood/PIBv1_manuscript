#!/bin/bash
#SBATCH -J PIBv1simsnf
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH -t 48:00:00
#SBATCH -e PIBv1_simulation_8_nf_runner_20251025_%j.stderr
#SBATCH -o PIBv1_simulation_8_nf_runner_20251025_%j.stdout
#SBATCH --mail-user patrick.f.reilly@yale.edu
#SBATCH --mail-type ALL
#SBATCH -p week

#Load the appropriate modules
module load nextflow/24.10.2b5932

pushd /gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes/PIBv1_ReviewerResponses/Simulations/Sim8/

#Name for the assembly run:
run_name="PIBv1"
#Name of the .config file in the current working directory:
cfg="${run_name}_simulation_8.config"
#Path to the Nextflow pipeline script:
pipelinedir="/gpfs/gibbs/pi/tucci/pfr8/bin/PIBv1_manuscript/Analysis_Pipelines"
pipeline="${pipelinedir}/adaptive_introgression_simulations.nf"
#Path to the Nextflow working directory, which I treat as scratch:
scratchdir="/gpfs/gibbs/pi/tucci/pfr8/nextflow_scratch"
scratch="${scratchdir}/${run_name}_sims_8/"
#Nextflow profiles to use (comma-separated string):
profiles="mccleary,hs37d5"
#Date for logs:
logdate="20251025"
#Resume?:
resume=""
resumenote=""

mkdir -p ${scratch}

NXF_OPTS="-Xms500M -Xmx3500M" /usr/bin/time -v nextflow -c ${cfg} run -ansi-log false ${pipeline} -profile ${profiles} ${resume} -w ${scratch} 2> ${run_name}_simulation_8_nextflow_${resumenote}${logdate}.stderr > ${run_name}_simulation_8_nextflow_${resumenote}${logdate}.stdout

popd


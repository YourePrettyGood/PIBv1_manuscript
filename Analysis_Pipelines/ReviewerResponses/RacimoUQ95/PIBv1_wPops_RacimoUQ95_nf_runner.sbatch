#!/bin/bash
#SBATCH -J PIBv1wPUQ95nf
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH -t 48:00:00
#SBATCH -e PIBv1_wPops_RacimoUQ95_nf_runner_20241231_%j.stderr
#SBATCH -o PIBv1_wPops_RacimoUQ95_nf_runner_20241231_%j.stdout
#SBATCH --mail-user patrick.f.reilly@yale.edu
#SBATCH --mail-type ALL
#SBATCH -p week

#Load the appropriate modules
module load nextflow/24.10.2b5932

pushd /gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes/PIBv1_ReviewerResponses/RacimoUQ95/

#Name for the assembly run:
run_name="PIBv1_wPops"
#Name of the .config file in the current working directory:
cfg="${run_name}_RacimoUQ95.config"
#Path to the Nextflow pipeline script:
pipelinedir="/gpfs/gibbs/pi/tucci/pfr8/bin/PIB_analysis_nextflow"
pipeline="${pipelinedir}/RacimoUQ95.nf"
#Path to the Nextflow working directory, which I treat as scratch:
scratchdir="/gpfs/gibbs/pi/tucci/pfr8/nextflow_scratch"
scratch="${scratchdir}/RacimoUQ95_${run_name}/"
#Nextflow profiles to use (comma-separated string):
profiles="mccleary,hs37d5"
#Date for logs:
logdate="20241231"
#Resume?:
resume=""
resumenote=""

mkdir -p ${scratch}

NXF_OPTS="-Xms500M -Xmx3500M" /usr/bin/time -v nextflow -c ${cfg} run -ansi-log false ${pipeline} -profile ${profiles} ${resume} -w ${scratch} 2> RacimoUQ95_${run_name}_nextflow_${resumenote}${logdate}.stderr > RacimoUQ95_${run_name}_nextflow_${resumenote}${logdate}.stdout

popd


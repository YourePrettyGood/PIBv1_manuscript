#!/bin/bash
#SBATCH -J AdIntFCSmerge
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32000
#SBATCH -t 24:00:00
#SBATCH -e logs/PIBv1_AdInt_FCS_mergenearbywindows_20251120_%j.stderr
#SBATCH -o logs/PIBv1_AdInt_FCS_mergenearbywindows_20251120_%j.stdout
#SBATCH --mail-user patrick.f.reilly@yale.edu
#SBATCH --mail-type ALL
#SBATCH -p day

#For mccleary:
module load R/4.3.2-foss-2022b

#Array of analysis groups to iterate over and run bedtools closest:
AGs=("Ata" "Baining-Kagat" "Baining-Mali" "Bellona-Rennell" "Kove" "Lavongai-Mussau" "Malaita" "Mamusi" "Melamela" "Nailik-Notsi-Tigak" "Nakanai-Mangseng" "Nasioi" "Santa-Cruz" "Saposa" "Sepik-Goroka" "Tikopia" "Vella-Lavella")
#Path to the BED of coding genes from GENCODE release 38 lifted to hs37d5:
annotbed="/gpfs/gibbs/pi/tucci/dt637/hg19_scikit-allel/key_files/gencode.v38lift37.annotation.codingGENES_sorted.bed"
#Path to the .genome index of hs37d5 to speed up bedtools:
genome="/gpfs/gibbs/pi/tucci/pfr8/refs/1kGP/hs37d5/hs37d5.genome"

pushd /gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes/PIBv1_ReviewerResponses/Simulations/

/usr/bin/time -v ./PIBv1_AdInt_FCS_mergenearbywindows_20251120.R

module purge
module load bedtools/b891a0b

for AG in ${AGs[@]};
   do
   hitbed="PIBv1_OCN_scans/PIBv1_${AG}_AdInt_overlapping_hits_neglog10pFCSge2_merged10k.bed"
   hitbed_wgenes=${hitbed//.bed/_wGenes.bed}
   /usr/bin/time -v bedtools closest -a <(sort -k1,1V -k2,2n -k3,3n ${hitbed}) -b ${annotbed} -g ${genome} -d -t all | \
      /usr/bin/time -v ./AdIntGeneList.awk > ${hitbed_wgenes}
done

popd

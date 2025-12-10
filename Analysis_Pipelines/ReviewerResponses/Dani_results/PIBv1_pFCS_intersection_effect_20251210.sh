#!/bin/bash

#This script does some simple reformatting of Dani's adaptive introgression
# pipeline outputs in preparation for generating histograms of p_FCS before
# and after overlapping with high frequency core haplotypes.

#Path to the high frequency core haplotypes with adjusted population
# labels:
highfreqcorehaps="/gpfs/gibbs/pi/tucci/dt637/hg19_scikit-allel/key_files/PIBv1_Sprime_putative_adaptive_introgressed_corehaps_OCNnoVanuatu_wGeneList_lumpedpopsrenamed.bed"
#Path to the header line for the output of the overlapping step:
outheader="/gpfs/gibbs/pi/tucci/dt637/hg19_scikit-allel/key_files/header_file.txt"

#Load the BEDtools module for the overlap step:
module load bedtools/b891a0b

#Symlink the relevant results files over:
ln -s /gpfs/gibbs/pi/tucci/dt637/hg19_scikit-allel/D_FCS_v3/*FCS.tsv .

#Do the overlap with high frequency core haplotypes from the same
# population:
#This step is effectively identical to what Dani does in
# run_bedtools_intersect_introgressed_tracks.sh, except
# done for the full *FCS.tsv files instead of the *FCS_0.01.tsv
# files (which are filtered).
#It takes about 1.5 minutes per population, maybe 2 at worst,
# and shouldn't take very much memory, though I allocate 8 GB just in case.
date
for i in *FCS.tsv;
   do
   fn=$(basename ${i});
   outfn=${fn//.tsv/_introgressed_corehaps_uniquePOP.tsv};
   echo "${outfn}";
   bedtools intersect -a <(tail -n+2 ${i}) -b <(tail -n+2 ${highfreqcorehaps}) -wao | \
      awk 'BEGIN{FS="\t";OFS=FS;}FNR==NR&&FNR==1{print;}FNR<NR&&$6==$22{print;}' ${outheader} - > ${outfn};
   date;
done

module unload bedtools/b891a0b

#Extract relevant columns and reformat selection scan windows for all
# populations:
#This step takes about 6 minutes per population, so is the longest step
# overall, but takes minimal memory (< 1 GB).
date
for i in *_FCS.tsv;
   do
   pop=${i%%_*};
   outfn="${pop}_PIBv1_selection_pvalues.tsv.gz";
   echo "${outfn}";
   ./reformatDTMresults.awk -v "overlapping=0" ${i} | gzip -9 > ${outfn};
   date;
done

#Extract relevant columns and reformat overlapping selection scan
# windows for all populations:
#This step takes 3-5 seconds per population and minimal memory (< 1 GB).
date
for i in *_uniquePOP.tsv;
   do
   pop=${i%%_*};
   outfn="${pop}_PIBv1_selection_highfreqcorehapoverlap_pvalues.tsv.gz";
   echo "${outfn}";
   ./reformatDTMresults.awk -v "overlapping=1" ${i} | gzip -9 > ${outfn};
   date;
done

#Now run the R script to generate the plot:
#Note that this uses ~29 GB RAM, but only takes ~6 minutes
# given 2 cores.
#The loading of input data takes a while (slight speedup
# if you allocate multiple cores, as readr's read_tsv()
# is multi-threaded), but so does the rendering of the
# plot, since the code is histogramming tens of millions
# of points per population.
module load R/4.3.2-foss-2022b
./PIBv1_pFCS_intersection_effect_20251210.R

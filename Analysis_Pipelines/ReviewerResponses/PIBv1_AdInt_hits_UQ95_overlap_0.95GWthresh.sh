#!/bin/bash

module load bedtools/b891a0b
UQ95thresh="AI_scan_hits/PIBv1_UQ95_0.95_thresholds.tsv"
genome="[path to ref]/hs37d5.genome"
overlap="AI_scan_hits/PIBv1_AdInt_scan_hit_UQ95_overlap_0.95GWthresh.tsv"
UQ95hits="RacimoUQ95/UQ95/PIBv1_wPops_UQ95_w40000.tsv"
UQ95_A=("AFR_YOR" "AFR_BNK,AFR_BNS,AFR_BIA,AFR_MAN,AFR_MBU,AFR_SAN,AFR_YOR" "AFR_.,EAS_." "AFR_.,EUR_.,EAS_.")

printf "TargetPop\tSourceMethod\tOverlap\tAvgHits\tNumHits\tNumCovered\tNumTracts\tOutgroup\n" > ${overlap}
echo "Header done"

for og in "${UQ95_A[@]}";
   do
   while read -a p;
      do
      AIhits="Simulations/PIBv1_OCN_scans/PIBv1_${p[0]}_AdInt_overlapping_hits_neglog10pFCSge2_merged10k.bed";
      printf "${p[0]}\tPIBv1\t" >> ${overlap};
      bedtools intersect -a <(cut -f1-4 ${AIhits}) -b <(./selectUQ95byBandThresh.awk -v "A=${og}" -v "B=${p[1]}" ${UQ95thresh} ${UQ95hits} | sort -k1,1V -k2,2n -k3,3n) -g ${genome} -sorted -wao | \
         ./tractCoveredSummary.awk -v "outgroup=${og}" >> ${overlap};
      echo "${p[0]} PIBv1 source ${og} outgroup done";
      printf "${p[0]}\tUQ95\t" >> ${overlap};
      bedtools intersect -a <(./selectUQ95byBandThresh.awk -v "A=${og}" -v "B=${p[1]}" ${UQ95thresh} ${UQ95hits} | sort -k1,1V -k2,2n -k3,3n) -b <(cut -f1-4 ${AIhits}) -g ${genome} -sorted -wao | \
         ./tractCoveredSummary.awk -v "outgroup=${og}" >> ${overlap};
      echo "${p[0]} UQ95 source ${og} outgroup done";
   done < AI_scan_hits/PIBv1_AI_scan_pop_map.tsv;
done

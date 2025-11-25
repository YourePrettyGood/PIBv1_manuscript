#!/bin/bash

module load bedtools/b891a0b

genome="[path to ref]/hs37d5.genome"
gaps="[path to ref]/hs37d5_assembly_gaps.bed"
deserts="deserts/PIBv1_archaic_deserts_hs37d5unmasked_maxgap0.bed"
introgressed="projections/PIBv1_total_tiling_paths_maxgap0.bed"
hist_out="Bstat/bStatistic_MurphyAndMcVicker_deserts_introgression_genomewide_hist_maxgap0.tsv"

Bstat_types=("McVicker" "Murphy_CADD" "Murphy_phastCons")

printf "Bstat\tcount\teCDF\tsource\tBstatType\n" > ${hist_out}
echo "Header done"

for t in "${Bstat_types[@]}";
   do
   bstat="Bstat/${t}_bStatistic_autosomes.bed"
   ./tractBstateCDF.awk -v "region=genomewide" -v "Btype=${t}" ${bstat} >> ${hist_out}
   echo "Genomewide ${t} done"
   bedtools intersect -a ${deserts} -b ${bstat} -g ${genome} -sorted -wao | \
      ./tractBstateCDF.awk -v "region=deserts" -v "Btype=${t}" >> ${hist_out}
   echo "Deserts ${t} done"
   bedtools complement -i ${deserts} -g ${genome} | \
      awk 'BEGIN{FS="\t";OFS=FS;}$1~/^[0-9][0-9]?$/{print;}' | \
      bedtools subtract -a - -b ${gaps} | \
      bedtools intersect -a - -b ${bstat} -g ${genome} -sorted -wao | \
      ./tractBstateCDF.awk -v "region=nondeserts" -v "Btype=${t}" >> ${hist_out}
   echo "Non-deserts ${t} done"
   bedtools intersect -a ${introgressed} -b ${bstat} -g ${genome} -sorted -wao | \
      ./tractBstateCDF.awk -v "region=introgressed" -v "Btype=${t}" >> ${hist_out}
   echo "Introgressed ${t} done"
   bedtools complement -i ${introgressed} -g ${genome} | \
      awk 'BEGIN{FS="\t";OFS=FS;}$1~/^[0-9][0-9]?$/{print;}' | \
      bedtools subtract -a - -b ${gaps} | \
      bedtools intersect -a - -b ${bstat} -g ${genome} -sorted -wao | \
      ./tractBstateCDF.awk -v "region=nonintrogressed" -v "Btype=${t}" >> ${hist_out}
   echo "Non-introgressed ${t} done"
done

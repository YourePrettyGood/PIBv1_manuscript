#!/bin/awk -f
#This script takes in the output of:
# bedtools intersect -a [per-indiv ROH tract BED] -b [per-indiv DP BEDGRAPH]
#  -wao -g [.genome file] -sorted
# and outputs the average DP along each ROH tract.
BEGIN{
   FS="\t";
   OFS=FS;
}
#Columns 1-3 define the ROH tract, 8 is the DP, 9 is the overlap length,
# so the average DP along the ROH is just sum_i 8_i*9_i / sum_i 9_i
{
   num[$1,$2,$3,$4]+=$8*$9;
   denom[$1,$2,$3,$4]+=$9;
}
END{
   #Output the average DP per ROH in some order, we can sort later:
   for (i in denom) {
      split(i, a, SUBSEP);
      printf "%s\t%i\t%i\t%s", a[1], a[2], a[3], a[4];
      if (denom[i] > 0) {
         printf "\t%f", num[i]/denom[i];
      } else {
         printf "\tNA";
      };
      printf "\n";
   };
}

#!/bin/awk -f
#This script takes the output of bedtools intersect -wao between
# A, a set of target tracts, and B, a set of query tracts, and
# calculates the proportion of target tracts covered by at least
# one query tract, the average number of query tracts covering
# a target tract, and the relevant numerators and denominator.
#The original use of this was to calculate how many adaptive
# introgression hits from the PIBv1 manuscript were reproduced
# by the Racimo et al. 2017 U and Q95 statistics.
#
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(outgroup) == 0) {
      print "Missing outgroup variable, please set it." > "/dev/stderr";
      exit 2;
   };
   #A file must have at least 3 columns, so default Bstart is 4:
   Bstart=4;
}
{
   #Auto-detect the start of the B file in the -wao output:
   if ($1 != $Bstart) {
      for (i=Bstart+1; i<=NF; i++) {
         if ($1 == $i || $i == ".") {
            Bstart=i;
            break;
         };
      };
   };
   if ($Bstart != ".") {
      #If the A tract has a hit in B, increment the hit count and
      # set the covered flag for this A tract:
      AIhit[$1,$2,$3]+=1;
      AIcovered[$1,$2,$3]=1;
   } else {
      #If the A tract doesn't have a hit in B, zero out the count
      # and covered flag for this A tract so it still is stored:
      AIhit[$1,$2,$3]=0;
      AIcovered[$1,$2,$3]=0;
   };
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   #Sum up the hits and covered flags and keep tract of the total
   # number of A tracts as denominator:
   hits=0;
   covered=0;
   tracts=0;
   for (t in AIhit) {
      hits+=AIhit[t];
      covered+=AIcovered[t];
      tracts+=1;
   };
   #Print the proportion of A tracts covered, average # hits per
   # A tract, and the numerators and denominator:
   if (tracts > 0) {
      print covered/tracts, hits/tracts, hits, covered, tracts, outgroup;
   } else {
      print "NA", "NA", hits, covered, tracts, outgroup;
   };
}

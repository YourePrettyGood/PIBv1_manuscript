#!/bin/awk -f
#This script does some simple reformatting of the output
# files from Dani's adaptive introgression pipeline. Two
# main input file formats are handled:
#1) *_FCS.tsv, which is the file of all selection scan windows with p_FCS
#   calculated for each (the default format, or overlapping <= 0)
#2) *_uniquePOP.tsv, which is the file of selection scan windows overlapping
#   high frequency core haplotypes from the same population (overlapping > 0)
#Differentiate between these two using the input argument "overlapping".
#I usually gzip the output of this to save some disk space.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(overlapping) == 0) {
      overlapping=0;
   };
}
NR==1{
   #Read in the input column names:
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
   #Output a header line:
   if (overlapping > 0) {
      print "chr", "start", "end", "population", "score_PBS", "p_PBS", "maxp_XPEHH", "p_HMP", "p_FCS", "CoreHapFreq";
   } else {
      print "chr", "start", "end", "population", "score_PBS", "p_PBS", "maxp_XPEHH", "p_HMP", "p_FCS";
   };
}
NR>1{
   #Handle NAs and convert blank values to NAs:
   if ($cols["score"] == "") {
      score="NA";
   } else {
      score=$cols["score"];
   };
   if ($cols["quantile_rank"] == "" || $cols["quantile_rank"] == "NA") {
      p_PBS="NA";
   } else {
      p_PBS=1-$cols["quantile_rank"];
   };
   if ($cols["max_quantile_rank_xpehh"] == "" || $cols["max_quantile_rank_xpehh"] == "NA") {
      p_XPEHH="NA";
   } else {
      p_XPEHH=1-$cols["max_quantile_rank_xpehh"];
   };
   #Output the window score, quantile ranks, and p_FCS, plus core hap
   # freq if provided:
   if (overlapping > 0) {
      print $cols["chr"], $cols["start"], $cols["end"], $cols["population"], score, p_PBS, p_XPEHH, $cols["HMP"], $cols["FCS"], $cols["TractFreq"];
   } else {
      print $cols["chr"], $cols["start"], $cols["end"], $cols["population"], score, p_PBS, p_XPEHH, $cols["HMP"], $cols["FCS"];
   };
}

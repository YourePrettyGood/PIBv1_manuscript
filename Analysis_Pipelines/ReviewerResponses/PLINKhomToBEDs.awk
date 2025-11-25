#!/bin/awk -f
#This script takes a PLINK .hom file from ROH analysis as input and
# outputs a bunch of BED files, one per sample found in the .hom file,
# containing the ROH tracts for that sample. This can then be used
# to calculate the average DP along each ROH tract.
BEGIN{
   OFS="\t";
   if (length(outprefix) == 0) {
      print "Missing outprefix, please provide it. Quitting." > "/dev/stderr";
      exit 2;
   };
}
NR==1{
   for (i=1; i<=NF; i++) {
      cols[$i]=i;
   };
}
NR>1{
   print $cols["CHR"], $cols["POS1"]-1, $cols["POS2"], $cols["IID"] > outprefix"_"$cols["IID"]"_ROH.bed";
}

#!/bin/awk -f
#This script generates an empirical cumulative density function for
# an input BED-like file. If the input is the output of bedtools intersect
# -wao, pass something other than "genomewide" to region, and this
# script will generate the eCDF specifically for intervals overlapping
# between the two inputs of the bedtools intersect call.
#This is achieved by using the overlap count in the last column
# and screening out records of -a that don't have overlap with -b.
#The original purpose of this script was to generate eCDFs for various
# versions of the B statistic of genomic constraint (e.g. the original
# B statistic from McVean, as well as the revised B statistics of
# Murphy et al. using CADD and phastCons) genome-wide as compared to
# only in archaic introgression deserts or archaic introgressed regions.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(region) == 0) {
      print "Missing region variable, please set this to the label for the tracts." > "/dev/stderr";
      exit 2;
   };
   if (length(Btype) == 0) {
      print "Missing Btype variable, please set this to the label for the B statistic source." > "/dev/stderr";
      exit 3;
   };
   #B statistic in principle goes from 0 to 1000:
   if (length(minval) == 0) {
      minval=0;
   };
   if (length(maxval) == 0) {
      maxval=1000;
   };
}
{
   if (region == "genomewide") {
      stat=$4;
      count=$3-$2;
   } else {
      stat=$(NF-1);
      count=$NF;
   };
   if (stat == ".") {
      next;
   };
   if (stat < minval || stat > maxval) {
      print "Value of statistic ("stat") appears to be outside the range "minval" <= B <= "maxval", please set minval and maxval correctly." > "/dev/stderr";
      exit 4;
   };
   hist[stat]+=count;
   sum+=count;
}
END{
   PROCINFO["sorted_in"]="@ind_num_asc";
   cumsum=0;
   if (complete > 0) {
      for (i=minval; i<=maxval; i++) {
         if (i in hist) {
            cumsum+=hist[i];
            print i, hist[i], cumsum/sum, region, Btype;
         } else {
            print i, 0, cumsum/sum, region, Btype;
         };
      };
   } else {
      for (b in hist) {
         cumsum+=hist[b];
         print b, hist[b], cumsum/sum, region, Btype;
      };
   };
}

#!/usr/bin/env nextflow
/* Pipeline to identify and characterize adaptive introgressed tracts from  *
 * archaic hominins (e.g. Neandertals and Denisovans)                       *
 * Core steps:                                                              *
 *  Extract archaic-matching S' variant sites and cluster by r^2 ->         *
 *  (Estimate tract frequencies && Find genes overlapping tracts)           */

//Default paths, globs, and regexes:
//Jointly genotyped VCFs:
params.vcf_glob = "${projectDir}/perchrom_VCFs/*.vcf.gz"
//Regex for extracting chromosome from VCF filename:
params.vcf_regex = ~/_chr(\p{Alnum}+)$/

//Sample ID file for exclusion:
params.samples_to_exclude = "<(echo)"

//Tract origins to assess and criteria to use:
params.origin_criteria = "PFR"
params.tract_origin_list = "Neandertal,Denisovan,Ambiguous"

//Reference-related parameters for the pipeline:
params.ref_prefix = "/gpfs/gibbs/pi/tucci/pfr8/refs"
params.ref = "${params.ref_prefix}/1kGP/hs37d5/hs37d5.fa"
params.autosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
params.nuclearchroms = "${params.autosomes},X,Y"

//Set up the channels of per-chromosome jointly-genotyped VCFs and their indices:
Channel
   .fromPath(params.vcf_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find VCFs matching glob: ${params.vcf_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from VCF: ${params.vcf_regex}" }
   .tap { perchrom_vcfs }
   .subscribe { println "Added ${it[1]} to perchrom_vcfs channel" }

Channel
   .fromPath(params.vcf_glob+'.tbi', checkIfExists: true)
   .ifEmpty { error "Unable to find VCF indices matching glob: ${params.vcf_glob}.tbi" }
   .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from VCF index: ${params.vcf_regex}" }
   .tap { perchrom_tbis }
   .subscribe { println "Added ${it[1]} to perchrom_tbis channel" }

//Set up the channel of target populations for Sprime:
Channel
   .fromPath(params.target_pops_file, checkIfExists: true)
   .ifEmpty { error "Unable to find target populations file for Sprime: ${params.target_pop_file}" }
   .splitCsv(sep: "\t")
   .map { it[0] }
   .tap { sprime_pops }
   .tap { sprime_target_pops }
   .tap { sprime_pops_for_tractfreqs }
   .subscribe { println "Added ${it} to sprime_pops channel" }

//Set up the channels of Sprime match and match rate files:
Channel
   .fromPath(params.sprime_match_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find Sprime match files matching glob: ${params.sprime_match_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.sprime_match_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract population ID from Sprime matches file: ${params.sprime_match_regex}" }
   .tap { sprime_matches }
   .tap { sprime_matches_for_score }
   .subscribe { println "Added ${it[1]} to sprime_matches channel" }
Channel
   .fromPath(params.sprime_matchrate_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find Sprime matchrate files matching glob: ${params.sprime_matchrate_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.sprime_matchrate_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract population ID from Sprime matchrate file: ${params.sprime_matchrate_regex}" }
   .tap { sprime_matchrates }
   .subscribe { println "Added ${it[1]} to sprime_matchrates channel" }

//Set up the channel of tract origins:
Channel
   .fromList(params.tract_origin_list.tokenize(','))
   .tap { tract_origins }
   .subscribe { println "Added ${it} to tract_origins channel" }

//Set up the file channels for the metadata files and the annotation GFF:
metadata = file(params.metadata_file, checkIfExists: true)
annotation_gff = file(params.annotation_gff, checkIfExists: true)
//And a file channel for the sample IDs to exclude:
excluded_samples = file(params.samples_to_exclude, checkIfExists: true)

//Default parameter values:
//Sample ID column name in metadata file:
params.id_colname = "SampleID"
//Column name for larger region encompassing target group for r^2 calcs:
params.r2_colname = "Region"
//Target group column name:
params.sprime_target_colname = "AnalysisGroup"
//Whether r^2 should be calculated based on genotypes or haplotypes:
params.phased = false
if (params.phased) {
   params.r2_type = "hap"
} else {
   params.r2_type = "geno"
}
//Minimum r^2 to report:
params.min_r2 = "0.3"
//Whether to run any of the projection code on haplotypes:
params.happroject = false

//Defaults for cpus, memory, and time for each process:
//Identify Sprime "tag" sites
params.tagsites_cpus = 1
params.tagsites_mem = 1
params.tagsites_timeout = '3h'
//Identify core haplotypes based on LD connected components
params.calcld_cpus = 1
params.calcld_mem = 4
params.calcld_timeout = '6h'
//Concatenate core haplotypes across chromosomes and tract origins
params.catcorehaps_cpus = 1
params.catcorehaps_mem = 1
params.catcorehaps_timeout = '3h'
//Annotate .score/_matches.tsv file with core haplotypes
params.corescore_cpus = 1
params.corescore_mem = 4
params.corescore_timeout = '3h'
//Sprime core haplotype projection
params.coreproject_cpus = 1
params.coreproject_mem = 1
params.coreproject_timeout = '24h'
//Sprime core haplotype frequencies
params.coretf_cpus = 1
params.coretf_mem = 1
params.coretf_timeout = '24h'
//Sprime core haplotype gene lists
params.coregenes_cpus = 1
params.coregenes_mem = 5
params.coregenes_timeout = '24h'
//Combine outputs across targets and chromosomes
params.catouts_cpus = 1
params.catouts_mem = 1
params.catouts_timeout = '2h'

//Preprocess the per-chromosome VCF channel to include the indices:
perchrom_vcfs
   .join(perchrom_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
   .filter({ params.autosomes.tokenize(',').contains(it[0]) })
   .tap { perchrom_vcfs_for_project }
   .tap { perchrom_vcfs_for_tractfreqs }
   .set { perchrom_vcfs_r2 }

process tagsites {
   tag "${pop} ${origin}"

   cpus params.tagsites_cpus
   memory { params.tagsites_mem.plus(task.attempt.minus(1).multiply(8))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.tagsites_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(pop), path(matches), path(matchrates), val(origin) from sprime_matches
      .join(sprime_matchrates, by: 0, failOnMismatch: true, failOnDuplicate: true)
      .combine(tract_origins)

   output:
   tuple val(pop), val(origin), path("${pop}_${params.origin_criteria}_${origin}Match_sites.tsv") into tagsites_output

   shell:
   '''
   #Identify archaic-matching S' sites based on tract match rate and
   # archaic matching:
   !{projectDir}/HumanPopGenScripts/Sprime/extract_Sprime_arcmatch_sites.awk -v "source=!{origin}" -v "criteria=!{params.origin_criteria}" -v "only_matches=1" <(gzip -dc !{matchrates}) <(gzip -dc !{matches}) > !{pop}_!{params.origin_criteria}_!{origin}Match_sites.tsv
   '''
}

process calc_ld {
   tag "${pop} ${origin} ${chrom}"

   cpus params.calcld_cpus
   memory { params.calcld_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.calcld_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/corehaps_logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), val(origin), path(sites), val(chrom), path(vcf), path(tbi) from tagsites_output
      .combine(perchrom_vcfs_r2)
   path metadata
   path excluded_samples

   output:
   tuple path("bcftools_view_${pop}_chr${chrom}_${origin}.stderr"), path("vcftools_${pop}_chr${chrom}_${params.r2_type}ld_r${params.min_r2}_${origin}.stderr"), path("vcftools_${pop}_chr${chrom}_${params.r2_type}ld_r${params.min_r2}_${origin}.stdout") into calc_ld_logs
   tuple val(pop), path("${pop}_chr${chrom}_Sprime_${origin}_${params.r2_type}ld_r${params.min_r2}_core_haplotypes.tsv.gz") into calc_ld_output

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_vcftools}
   module load !{params.mod_R}
   #Identify the set of samples from the superpopulation containing the
   # target population, since we want LD calculated from a large sample:
   superpop=$(!{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.r2_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{pop}" !{metadata} | uniq)
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.r2_colname}" -v "select=${superpop}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_superpop_samples.tsv
   #Pull out the archaic-matching S' sites and calculate pairwise r^2:
   awk -v "chrom=!{chrom}" '\$1==chrom' !{sites} > !{origin}_sites.tsv
   if [[ -s "!{pop}_superpop_samples.tsv" ]] && [[ -s "!{origin}_sites.tsv" ]]; then
      bcftools view -R !{origin}_sites.tsv -T !{origin}_sites.tsv -S !{pop}_superpop_samples.tsv -a -Ov !{vcf} 2> bcftools_view_!{pop}_chr!{chrom}_!{origin}.stderr | \
         vcftools --vcf - --out !{pop}_chr!{chrom}_!{origin}_LD_trimalt --!{params.r2_type}-r2 --min-r2 !{params.min_r2} 2> vcftools_!{pop}_chr!{chrom}_!{params.r2_type}ld_r!{params.min_r2}_!{origin}.stderr > vcftools_!{pop}_chr!{chrom}_!{params.r2_type}ld_r!{params.min_r2}_!{origin}.stdout
      #Filter and reformat the r^2 results to only between pairs of input sites:
      !{projectDir}/HumanPopGenScripts/Sprime/extract_arcmatch_tag_SNPs_Rsquared.awk -v "source=vcftools" !{origin}_sites.tsv !{pop}_chr!{chrom}_!{origin}_LD_trimalt.!{params.r2_type}.ld > !{pop}_chr!{chrom}_!{origin}_!{params.r2_type}ld_r!{params.min_r2}_rsquared.tsv
      #Construct a graph from these r^2 edges and identify core haplotypes as
      # connected components of the graph, outputting a labeled list of sites:
      !{projectDir}/HumanPopGenScripts/Sprime/LD_connected_components.R !{pop}_chr!{chrom}_!{origin}_!{params.r2_type}ld_r!{params.min_r2}_rsquared.tsv !{pop}_chr!{chrom}_Sprime_!{origin}_!{params.r2_type}ld_r!{params.min_r2}_core_haplotypes.tsv.gz !{chrom}
   else
     #Catch the case where no sites or samples are found and just pass through
     # a blank core haplotypes file plus logs indicating skipping the above steps:
     printf "" | \
        gzip -9 > !{pop}_chr!{chrom}_Sprime_!{origin}_!{params.r2_type}ld_r!{params.min_r2}_core_haplotypes.tsv.gz
     echo "Not enough samples or sites for !{pop} chr!{chrom} !{origin} (superpop=${superpop}), skipping LD clustering of tag SNPs." | \
        tee bcftools_view_!{pop}_chr!{chrom}_!{origin}.stderr | \
        tee vcftools_!{pop}_chr!{chrom}_!{params.r2_type}ld_r!{params.min_r2}_!{origin}.stderr | \
        tee vcftools_!{pop}_chr!{chrom}_!{params.r2_type}ld_r!{params.min_r2}_!{origin}.stdout
   fi
   '''
}

process cat_corehaps {
   tag "${pop}"

   cpus params.catcorehaps_cpus
   memory { params.catcorehaps_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.catcorehaps_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/corehaps_logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), path(corehaps) from calc_ld_output
      .collectFile() { [ "${it[0]}.txt", it[1].toAbsolutePath().toString()+'\t'+it[1].getName()+'\n' ] }
      .map({ [ it.getSimpleName(), it ] })

   output:
   tuple val(pop), path("${pop}_Sprime_${params.r2_type}ld_r${params.min_r2}_core_haplotypes_annotated.tsv.gz") into cat_corehaps_output

   shell:
   '''
   #Stage the core haplotype files:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[0]} ${a[1]};
   done < !{corehaps}
   #Concatenate the core haplotypes across both chromosomes and tract origins:
   cut -f2 !{corehaps} | \
      sort -t"_" -k1,1 -k2,2V -k4,4 | \
      while read a;
         do
         fn=$(basename ${a});
         prefix=${fn%_!{params.r2_type}ld_r!{params.min_r2}_core_haplotypes.tsv.gz};
         origin=${prefix#!{pop}_chr*_Sprime_};
         gzip -dc ${a} | \
            !{projectDir}/HumanPopGenScripts/Sprime/cat_core_haplotypes.awk -v "pop=!{pop}" -v "origin=${origin}";
      done | \
         gzip -9 > !{pop}_Sprime_!{params.r2_type}ld_r!{params.min_r2}_core_haplotypes_annotated.tsv.gz
   '''
}

process core_score {
   tag "${pop}"

   cpus params.corescore_cpus
   memory { params.corescore_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.corescore_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(pop), path(corehaps), path(matches) from cat_corehaps_output
      .join(sprime_matches_for_score, by: 0, failOnDuplicate: true, failOnMismatch: true)

   output:
   tuple val(pop), path("${pop}_Sprime_corehaps_matches.tsv") into corehaps_score_for_project,corehaps_score_for_tractfreqs

   shell:
   '''
   gzip -dc !{corehaps} | \
      !{projectDir}/HumanPopGenScripts/Sprime/applyCoreHaplotypesToScore.awk - <(gzip -dc !{matches}) > !{pop}_Sprime_corehaps_matches.tsv
   '''
}

process corehaps_project {
   tag "${pop} chr${chrom}"

   cpus params.coreproject_cpus
   memory { params.coreproject_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.coreproject_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/corehaps", mode: 'copy', pattern: '*.bed'

   input:
   tuple val(pop), path(corehapscores), val(chrom), path(vcf), path(tbi) from corehaps_score_for_project
      .combine(perchrom_vcfs_for_project)
   path metadata
   path excluded_samples

   output:
   tuple val(pop), path("${pop}_chr${chrom}_Sprime_corehaps_perSample.bed") into corehaps_project_BEDs

   shell:
   '''
   module load !{params.mod_bcftools}
   #Identify the samples for this population:
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{pop}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_pop_samples.tsv
   #Project Sprime alleles onto each individual to identify tracts from genotypes:
   bcftools query -f '%CHROM:%POS[\t%GT]\n' -r !{chrom} -H -S !{pop}_pop_samples.tsv !{vcf} | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimePerSampleTracts.awk !{corehapscores} - | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBED.awk | \
         sort -k1,1V -k2,2n -k3,3n > !{pop}_chr!{chrom}_Sprime_corehaps_perSample.bed
   '''
}

process corehaps_tractfreqs {
   tag "${pop} ${qpop} chr${chrom}"

   cpus params.coretf_cpus
   memory { params.coretf_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.coretf_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(qpop), val(pop), path(corehapscores), val(chrom), path(vcf), path(tbi) from sprime_pops_for_tractfreqs
      .combine(corehaps_score_for_tractfreqs)
      .combine(perchrom_vcfs_for_tractfreqs)
   path metadata
   path excluded_samples

   output:
   tuple val(pop), val(qpop), val(chrom), path("${pop}_chr${chrom}_Sprime_${qpop}_corehap_freqs.tsv.gz") into corehaps_tract_freqs_for_genes,corehaps_tract_freqs_for_cat

   shell:
   '''
   module load !{params.mod_bcftools}
   #Identify the samples for this population:
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{qpop}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{qpop}_pop_samples.tsv
   #Calculate the median Sprime allele frequency for each tract:
   bcftools query -f '%CHROM:%POS[\t%GT]\n' -r !{chrom} -H -S !{qpop}_pop_samples.tsv !{vcf} | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeArchaicAF.awk -v "spop=!{pop}" -v "pop=!{qpop}" -v "all=1" !{corehapscores} - | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractMedianAF.awk | \
      gzip -9 > !{pop}_chr!{chrom}_Sprime_!{qpop}_corehap_freqs.tsv.gz
   '''
}

process corehaps_genes {
   tag "${pop} chr${chrom}"

   cpus params.coregenes_cpus
   memory { params.coregenes_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.coregenes_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(pop), val(qpop), val(chrom), path("${pop}_chr${chrom}_Sprime_${qpop}_corehap_freqs.tsv.gz") from corehaps_tract_freqs_for_genes
      .filter({ it[0] == it[1] })
   path annotation_gff

   output:
   tuple val(pop), path("${pop}_chr${chrom}_Sprime_corehaps_gene_lists.tcsv.gz"), path("${pop}_chr${chrom}_Sprime_corehaps_gene_name_lists.tcsv.gz") into corehaps_gene_lists

   shell:
   '''
   module load !{params.mod_bedtools}
   #Identify genes overlapping Sprime tracts:
   bedtools intersect \
      -a <(gzip -dc !{pop}_chr!{chrom}_Sprime_!{pop}_corehap_freqs.tsv.gz | \
         awk 'BEGIN{FS="\t";OFS=FS;}{print $1, $2, $3, "TractID="$4";AAF="$6, ".", "+";}' | \
         sort -k1,1 -k2,2n -k3,3n) \
      -b <(gzip -dc !{annotation_gff} | \
         awk -v "keepchr=0" -v "chrom=!{chrom}" 'BEGIN{FS="\t";OFS=FS;if (length(keepchr) == 0) {keepchr=0;};}/^#/{print;}!/^#/{if (!keepchr) {sub("chr", "", $1);}; if ($1 == chrom) {print $0;};}') \
      -wao | \
      tee >(!{projectDir}/HumanPopGenScripts/Sprime/SprimeTractGeneList.awk -v "trim=1" -v "tag=gene_name" | \
         sort -k1,1V -k2,2n -k3,3n | \
         gzip -9 > !{pop}_chr!{chrom}_Sprime_corehaps_gene_name_lists.tcsv.gz) | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractGeneList.awk -v "trim=1" -v "tag=ID" | \
      sort -k1,1V -k2,2n -k3,3n | \
      gzip -9 > !{pop}_chr!{chrom}_Sprime_corehaps_gene_lists.tcsv.gz
   '''
}

process cat_outs {
   cpus params.catouts_cpus
   memory { params.catouts_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.catouts_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/corehaps", mode: 'copy', pattern: '*.t{sv,csv}.gz'

   input:
   path projections from corehaps_project_BEDs
      .map({ it[1] })
      .collectFile() { [ "projection_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path tractfreqs from corehaps_tract_freqs_for_cat
      .map({ it[3] })
      .collectFile() { [ "tractfreq_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path genelists from corehaps_gene_lists
      .map({ [ it[1], it[2] ] })
      .collectFile() { [ "genelist_paths.tsv", it[0].getSimpleName()+'\t'+it[0].getName()+'\t'+it[0]+'\n'+it[1].getSimpleName()+'\t'+it[1].getName()+'\t'+it[1]+'\n' ] }
   path targetpops from sprime_target_pops
      .collectFile(name: 'Sprime_target_populations.txt', newLine: true, sort: true)
   path metadata
   path excluded_samples

   output:
   tuple path("${params.run_name}_Sprime_perChrom_perIndiv_perPop_corehap_lengths.tsv.gz"), path("${params.run_name}_Sprime_allPopPairs_corehap_freqs.tsv.gz"), path("${params.run_name}_Sprime_targetpop_corehap_freqs.tsv.gz"), path("${params.run_name}_Sprime_corehaps_gene_name_lists.tcsv.gz"), path("${params.run_name}_Sprime_corehaps_gene_lists.tcsv.gz") into catouts_outputs

   shell:
   '''
   #Set the array of chromosomes to iterate over:
   IFS="," read -a chroms <<< "!{params.autosomes}"
   #Symlink the many input files to avoid SLURM SBATCH script size limits:
   cat projection_paths.tsv tractfreq_paths.tsv genelist_paths.tsv | \
      while IFS=$'\t' read -a a;
         do
         ln -s ${a[2]} ${a[1]};
      done
   #Calculate total per-individual tract length and
   # concatenate files across populations:
   header=1
   while IFS=$'\t' read p;
      do
      !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=${p}" !{metadata} | \
         !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBEDtoLengths.awk -v "header=${header}" -v "pop=${p}" - ${p}_chr{!{params.autosomes}}_Sprime_corehaps_perSample.bed;
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_perChrom_perIndiv_perPop_corehap_lengths.tsv.gz
   unset header
   #Concatenate tract frequency files and add a header:
   header=1
   for chrom in ${chroms[@]};
      do
      while IFS=$'\t' read p;
         do
         if [[ "${header}" == "1" ]]; then
            printf "Chromosome\tStart\tEnd\tTractID\tQueryPop\tMedianAF\tMinAF\tMaxAF\tNumSites\n"
            header=""
         fi
         while IFS=$'\t' read q;
            do
            gzip -dc ${p}_chr${chrom}_Sprime_${q}_corehap_freqs.tsv.gz | \
               sort -k1,1V -k2,2n -k3,3n;
         done < !{targetpops};
      done < !{targetpops};
   done | \
      gzip -9 > !{params.run_name}_Sprime_allPopPairs_corehap_freqs.tsv.gz
   unset header
   #Concatenate tract frequency files only for target pops as query and add a header:
   header=1
   while IFS=$'\t' read p;
      do
      for chrom in ${chroms[@]};
         do
         if [[ "${header}" == "1" ]]; then
            printf "Chromosome\tStart\tEnd\tTractID\tQueryPop\tMedianAF\tMinAF\tMaxAF\tNumSites\n"
            header=""
         fi
         gzip -dc ${p}_chr${chrom}_Sprime_${p}_corehap_freqs.tsv.gz | \
            sort -k1,1V -k2,2n -k3,3n;
      done;
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_targetpop_corehap_freqs.tsv.gz
   unset header
   #Concatenate the gene name lists across populations:
   header=1
   while IFS=$'\t' read p;
      do
      if [[ "${header}" == "1" ]]; then
         printf "Chromosome\tStart\tEnd\tTractID\tOverlappingGenes\n"
         header=""
      fi
      gzip -dc ${p}_chr{!{params.autosomes}}_Sprime_corehaps_gene_name_lists.tcsv.gz;
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_corehaps_gene_name_lists.tcsv.gz
   unset header
   #Concatenate the gene ID lists across populations:
   header=1
   while IFS=$'\t' read p;
      do
      if [[ "${header}" == "1" ]]; then
         printf "Chromosome\tStart\tEnd\tTractID\tOverlappingGenes\n"
         header=""
      fi
      gzip -dc ${p}_chr{!{params.autosomes}}_Sprime_corehaps_gene_lists.tcsv.gz;
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_corehaps_gene_lists.tcsv.gz
   unset header
   '''
}

/*
if (params.happroject) {
   sprime_matches
      .tap { sprime_matches_for_happroject }
   perchrom_vcfs_for_project
      .tap { perchrom_vcfs_for_happroject }
} else {
   sprime_matches_for_happroject = Channel.empty()
   perchrom_vcfs_for_happroject = Channel.empty()
}
process hap_project {
   tag "${pop} chr${chrom}"

   cpus params.happroject_cpus
   memory { params.happroject_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.happroject_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/happroject", mode: 'copy', pattern: '*.bed'

   input:
   tuple val(pop), path(matches), val(chrom), path(vcf), path(tbi) from sprime_matches_for_happroject
      .combine(perchrom_vcfs_for_happroject)
   path metadata
   path excluded_samples

   output:
   tuple val(pop), path("${pop}_chr${chrom}_Sprime_tracts_perSample_haplotypes.bed") into haps_project_BEDs

   shell:
   '''
   module load !{params.mod_bcftools}
   #Identify the samples for this population:
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{pop}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_pop_samples.tsv
   #Project Sprime alleles onto each individual to identify tracts from genotypes:
   bcftools query -f '%CHROM:%POS[\t%GT]\n' -r !{chrom} -H -S !{pop}_pop_samples.tsv !{vcf} | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimePerSampleTracts.awk -v "phased=1" !{matches} - | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBED.awk -v "phased=1" | \
         sort -k1,1V -k2,2n -k3,3n > !{pop}_chr!{chrom}_Sprime_tracts_perSample_haplotypes.bed
   '''
}
*/

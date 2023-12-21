#!/usr/bin/env nextflow
/* Pipeline to identify and characterize introgression tracts from archaic  *
 *  hominins (e.g. Neandertals and Denisovans)                              *
 * Core steps:                                                              *
 *  Extraction of target populations + outgroup from main VCF for Sprime -> *
 *  Concatenation of per-autosome VCFs for Sprime ->                        *
 *  Sprime for each target population ->                                    *
 *  (Calculate match rates && Estimate tract frequencies &&                 *
 *   Project S' tracts onto individuals && Find genes overlapping tracts)   */

//Default paths, globs, and regexes:
//Jointly genotyped VCFs:
params.vcf_glob = "${projectDir}/perchrom_VCFs/*.vcf.gz"
//Regex for extracting chromosome from VCF filename:
params.vcf_regex = ~/_chr(\p{Alnum}+)$/
//Recombination rate maps:
params.recmap_glob = "${projectDir}/"
//Regex for extracting chromosome from recombination rate map filename:
params.recmap_regex = ~/chr([0-9XY_par]+)/
//Archaic VCFs:
params.arcvcf_glob = "${projectDir}/archaic_VCFs/*.vcf.gz"

//Include/filter expression to apply to input VCFs:
params.input_filter_str = ""
//Sample ID file for exclusion:
params.samples_to_exclude = "<(echo)"

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

//Set up the channel of recombination rate map files for Sprime:
Channel
   .fromPath(params.recmap_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find recombination rate map files matching glob: ${params.recmap_glob}" }
   .map { a -> [ (a =~ params.recmap_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from VCF: ${params.recmap_glob} =~ ${params.recmap_regex}" }
   .tap { recmaps }
   .subscribe { println "Added chr${it[0]} recombination rate map (${it[1]}) to recmaps channel" }

//Set up the channel of target populations for Sprime:
Channel
   .fromPath(params.target_pops_file, checkIfExists: true)
   .ifEmpty { error "Unable to find target populations file for Sprime: ${params.target_pop_file}" }
   .splitCsv(sep: "\t")
   .map { it[0] }
   .tap { sprime_pops }
   .tap { sprime_target_pops }
   .subscribe { println "Added ${it} to sprime_pops channel" }

//Set up the channels of archaic VCFs and their indices:
Channel
   .fromPath(params.arcvcf_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find archaic VCFs matching glob: ${params.arcvcf_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from archaic VCF: ${params.vcf_regex}" }
   .tap { archaic_vcfs }
   .subscribe { println "Added ${it[1]} to archaic_vcfs channel" }

Channel
   .fromPath(params.arcvcf_glob+'.tbi', checkIfExists: true)
   .ifEmpty { error "Unable to find archaic VCF indices matching glob: ${params.arcvcf_glob}.tbi" }
   .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from archaic VCF index: ${params.vcf_regex}" }
   .tap { archaic_tbis }
   .subscribe { println "Added ${it[1]} to archaic_tbis channel" }

num_autosomes = params.autosomes.tokenize(',').size()

//Set up the file channels for the metadata files and the annotation GFF:
metadata = file(params.metadata_file, checkIfExists: true)
archaic_metadata = file(params.arc_metadata_file, checkIfExists: true)
annotation_gff = file(params.annotation_gff, checkIfExists: true)
//And a file channel for the sample IDs to exclude:
excluded_samples = file(params.samples_to_exclude, checkIfExists: true)

//Default parameter values:
//Sample ID column name in metadata file:
params.id_colname = "Sample"
//Outgroup name to use:
params.sprime_outgroup_colname = "Region"
params.sprime_outgroup = "Africa"
//Target group column name:
params.sprime_target_colname = "AnalysisGroup"

//Defaults for cpus, memory, and time for each process:
//VCF subsetting for Sprime
params.sprimesubset_cpus = 1
params.sprimesubset_mem = 4
params.sprimesubset_timeout = '24h'
params.sprimesubset_queue = 'day'
params.sprimesubset_long_queue = 'day'
//VCF concatenating for Sprime
params.concatvcf_cpus = 1
params.concatvcf_mem = 4
params.concatvcf_timeout = '24h'
params.concatvcf_queue = 'day'
params.concatvcf_long_queue = 'week'
//Sprime
params.sprime_cpus = 1
params.sprime_mem = 8
params.sprime_timeout = '24h'
params.sprime_queue = 'day'
params.sprime_long_queue = 'day'
//Sprime match rates
params.sprimematch_cpus = 1
params.sprimematch_mem = 4
params.sprimematch_timeout = '24h'
params.sprimematch_queue = 'day'
params.sprimematch_long_queue = 'day'
//Sprime projection
params.sprimeproject_cpus = 1
params.sprimeproject_mem = 1
params.sprimeproject_timeout = '24h'
params.sprimeproject_long_timeout = '48h'
params.sprimeproject_queue = 'day'
params.sprimeproject_long_queue = 'week'
//Sprime tract frequencies
params.sprimetf_cpus = 1
params.sprimetf_mem = 1
params.sprimetf_timeout = '24h'
params.sprimetf_long_timeout = '72h'
params.sprimetf_queue = 'day'
params.sprimetf_long_queue = 'week'
//Sprime tract gene lists
params.sprimegenes_cpus = 1
params.sprimegenes_mem = 5
params.sprimegenes_timeout = '24h'
params.sprimegenes_queue = 'day'
params.sprimegenes_long_queue = 'day'
//Combine outputs across targets and chromosomes
params.catouts_cpus = 1
params.catouts_mem = 1
params.catouts_timeout = '2h'
params.catouts_queue = 'day'
params.catouts_long_queue = 'day'

//Detection mechanism for "chr" prefixes of autosomes:
has_chr_prefix = params.autosomes.count('chr') > 0
autosome_list = params.autosomes
autosome_num_list = has_chr_prefix ? params.autosomes.replaceAll('chr', '') : autosome_list
Map autosome_map = [params.autosomes.replaceAll('chr', '').tokenize(','), params.autosomes.tokenize(',')]
   .transpose()
   .collectEntries()

//Preprocess the per-chromosome VCF channel to include the indices:
perchrom_vcfs
   .join(perchrom_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
   .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
   .set { perchrom_vcfs_subset }

//Do the same for the archaic VCF channel:
archaic_vcfs
   .join(archaic_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
   .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
   .flatMap({ [it[1], it[2]] })
   .collect()
   .set { archaic_vcfs_autosomes }

//Retain only the autosomes for the recombination rate maps:
recmaps
   .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
   .set { recmaps_autosomes }

process sprime_vcf_subset {
   tag "${pop} chr${chrom}"

   cpus params.sprimesubset_cpus
   memory { params.sprimesubset_mem.plus(task.attempt.minus(1).multiply(8))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.sprimesubset_timeout }
   queue { task.attempt >= 2 ? params.sprimesubset_long_queue : params.sprimesubset_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*_samples.tsv', saveAs: { 'chr'+chrom+'_'+it }

   input:
   tuple val(pop), val(chrom), path(input_vcf), path(input_tbi) from sprime_pops.combine(perchrom_vcfs_subset)
//   each pop from sprime_pops
   path metadata
   path excluded_samples

   output:
   tuple path("bcftools_view_selectpops_Sprime_${pop}_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_Sprime_${pop}_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_Sprime_${pop}_chr${chrom}.stdout") into sprime_vcf_subset_logs
   tuple val(pop), path("${pop}_chr${chrom}_nomissinggenos.vcf.gz"), path("${pop}_chr${chrom}_nomissinggenos.vcf.gz.tbi") into Sprime_perchrom_vcfs,Sprime_perchrom_vcfs_project,Sprime_perchrom_vcfs_tractfreqs
   tuple val(pop), path("${pop}_outgroup_samples.tsv") into Sprime_outgroup,Sprime_outgroup_project
   tuple val(pop), path("${pop}_samples.tsv") into Sprime_pop_samples

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_outgroup_colname}" -v "select=!{params.sprime_outgroup}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_outgroup_samples.tsv
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{pop}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_samples.tsv
   cat !{pop}_samples.tsv !{pop}_outgroup_samples.tsv > !{pop}_plusOutgroup_samples.tsv
   #We have to separate the sample selection from the genotype missingness filter
   # because of a bug in bcftools view when -S and -g are used in tandem:
   #Way too many sites get filtered if you use them together, possibly more than
   # if -g was used alone.
   bcftools view -S !{pop}_plusOutgroup_samples.tsv -Ou !{input_vcf} 2> bcftools_view_selectpops_Sprime_!{pop}_chr!{chrom}.stderr | \
   bcftools view -i '!{params.input_filter_str}' -g ^miss -Oz -o !{pop}_chr!{chrom}_nomissinggenos.vcf.gz 2> bcftools_view_nomissinggenos_Sprime_!{pop}_chr!{chrom}.stderr > bcftools_view_nomissinggenos_Sprime_!{pop}_chr!{chrom}.stdout
   tabix -f !{pop}_chr!{chrom}_nomissinggenos.vcf.gz
   '''
}

process concat_input_vcf {
   tag "${pop}"

   cpus params.concatvcf_cpus
   memory { params.concatvcf_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.concatvcf_timeout }
   queue { task.attempt >= 2 ? params.concatvcf_long_queue : params.concatvcf_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), path("*"), path("*") from Sprime_perchrom_vcfs.groupTuple(by: 0, size: num_autosomes)

   output:
   tuple path("bcftools_concat_Sprime_${pop}.stderr"), path("bcftools_concat_Sprime_${pop}.stdout") into concat_input_vcf_logs
   tuple val(pop), path("${pop}_nomissinggenos.vcf.gz"), path("${pop}_nomissinggenos.vcf.gz.tbi") into Sprime_vcfs

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   bcftools concat --threads !{task.cpus} -Oz -o !{pop}_nomissinggenos.vcf.gz !{pop}_chr{!{autosome_num_list}}_nomissinggenos.vcf.gz 2> bcftools_concat_Sprime_!{pop}.stderr > bcftools_concat_Sprime_!{pop}.stdout
   tabix -f !{pop}_nomissinggenos.vcf.gz
   '''
}

process sprime {
   tag "${pop} chr${chrom}"

   cpus params.sprime_cpus
   memory { params.sprime_mem.plus(1).plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.sprime_timeout }
   queue { task.attempt >= 2 ? params.sprime_long_queue : params.sprime_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.log'
   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.score'

   input:
   tuple val(pop), path("${pop}_nomissinggenos.vcf.gz"), path("${pop}_nomissinggenos.vcf.gz.tbi"), path("${pop}_outgroup_samples.tsv"), val(chrom), path(recmap) from Sprime_vcfs
      .join(Sprime_outgroup
         .groupTuple(by: 0, size: num_autosomes)
         .map({ [it[0], it[1].take(1)] }), by: 0, failOnDuplicate: true, failOnMismatch: true)
      .combine(recmaps_autosomes)

   output:
   tuple path("Sprime_${pop}_chr${chrom}.stderr"), path("Sprime_${pop}_chr${chrom}.stdout"), path("${pop}_chr${chrom}_Sprime.log") into Sprime_logs
   tuple val(pop), path("${pop}_chr${chrom}_Sprime.score") into Sprime_scores_matchrate,Sprime_scores_project,Sprime_scores_tractfreqs

   shell:
   sprime_retry_mem = params.sprime_mem.plus(task.attempt.minus(1).multiply(16))
   chrid = autosome_map[chrom]
   '''
   module load !{params.mod_Sprime}
   java -Xmx!{sprime_retry_mem}g -jar ${SPRIME} gt=!{pop}_nomissinggenos.vcf.gz outgroup=!{pop}_outgroup_samples.tsv map=!{recmap} chrom=!{chrid} out=!{pop}_chr!{chrom}_Sprime 2> Sprime_!{pop}_chr!{chrom}.stderr > Sprime_!{pop}_chr!{chrom}.stdout
   '''
}

process sprime_matchrate {
   tag "${pop}"

   cpus params.sprimematch_cpus
   memory { params.sprimematch_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.sprimematch_timeout }
   queue { task.attempt >= 2 ? params.sprimematch_long_queue : params.sprimematch_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

//   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.tsv.gz'

   input:
   tuple val(pop), path("*") from Sprime_scores_matchrate.groupTuple(by: 0, size: num_autosomes)
   path("*") from archaic_vcfs_autosomes
   path archaic_metadata

   output:
   path("${pop}_autosomes_Sprime_matches.tsv.gz") into Sprime_matches
   path("${pop}_autosomes_Sprime_match_rates.tsv.gz") into Sprime_matchrates

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_bedtools}
   #Identify Sprime alleles that match or mismatch the archaics:
   for chr in {!{autosome_num_list}};
      do
      #Establish the set of Sprime sites for this population and chromosome:
      fgrep -v "CHROM" !{pop}_chr${chr}_Sprime.score | \
         cut -f1,2 > !{pop}_chr${chr}_Sprime_sites.tsv
      !{projectDir}/HumanPopGenScripts/Sprime/archaicMatchSprime.awk -v "spop=!{pop}" !{archaic_metadata} \
         <(bcftools view -R !{pop}_chr${chr}_Sprime_sites.tsv -T !{pop}_chr${chr}_Sprime_sites.tsv *_chr${chr}.vcf.gz) \
         !{pop}_chr${chr}_Sprime.score
   done | \
      awk 'BEGIN{FS="\t";OFS=FS;}/^CHROM/&&NR==1{print;}!/^CHROM/{print;}' | \
      gzip -9 > !{pop}_autosomes_Sprime_matches.tsv.gz
   #Calculate the match rates from this file:
   gzip -dc !{pop}_autosomes_Sprime_matches.tsv.gz | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeArchaicMatchRate.awk | \
      gzip -9 > !{pop}_autosomes_Sprime_match_rates.tsv.gz
   '''
}

process sprime_project {
   tag "${pop}"

   cpus params.sprimeproject_cpus
   memory { params.sprimeproject_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? params.sprimeproject_long_timeout : params.sprimeproject_timeout }
   queue { task.attempt >= 2 ? params.sprimeproject_long_queue : params.sprimeproject_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

//   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.bed'

   input:
   tuple val(pop), path("*"), path("*"), path("*"), path("*") from Sprime_scores_project
      .groupTuple(by: 0, size: num_autosomes)
      .join(Sprime_perchrom_vcfs_project
         .groupTuple(by: 0, size: num_autosomes), by: 0, failOnDuplicate: true, failOnMismatch: true)
      .join(Sprime_outgroup_project
         .groupTuple(by: 0, size: num_autosomes)
         .map({ [it[0], it[1].take(1)]}), by: 0, failOnDuplicate: false, failOnMismatch: true)

   output:
   path("${pop}_Sprime_tracts_perSample.bed") into Sprime_project_BEDs

   shell:
   '''
   module load !{params.mod_bcftools}
   #Project Sprime alleles onto each individual to identify tracts from genotypes:
   for chr in {!{autosome_list}};
      do
      if [[ "!{has_chr_prefix}" == "true" ]]; then
         chrnum=${chr#chr};
      else
         chrnum=${chr};
      fi;
      bcftools query -f '%CHROM:%POS[\t%GT]\n' -r ${chr} -H -S ^!{pop}_outgroup_samples.tsv !{pop}_chr${chrnum}_nomissinggenos.vcf.gz | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimePerSampleTracts.awk !{pop}_chr${chrnum}_Sprime.score - | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBED.awk | \
         sort -k1,1V -k2,2n -k3,3n
   done > !{pop}_Sprime_tracts_perSample.bed
   '''
}

process sprime_tractfreqs {
   tag "${pop} ${qpop}"

   cpus params.sprimetf_cpus
   memory { params.sprimetf_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? params.sprimetf_long_timeout : params.sprimetf_timeout }
   queue { task.attempt >= 2 ? params.sprimetf_long_queue : params.sprimetf_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

//   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.tsv.gz'

   input:
   tuple val(pop), path("*"), val(qpop), path("*"), path("*"), path("*") from Sprime_scores_tractfreqs
      .groupTuple(by: 0, size: num_autosomes)
      .combine(Sprime_pop_samples
         .groupTuple(by: 0, size: num_autosomes)
         .map({ [it[0], it[1].take(1)] })
         .join(Sprime_perchrom_vcfs_tractfreqs
      .groupTuple(by: 0, size: num_autosomes), by: 0, failOnDuplicate: true, failOnMismatch: true))

   output:
   tuple val(pop), val(qpop), path("${pop}_Sprime_${qpop}_tract_freqs.tsv.gz") into Sprime_tract_freqs,Sprime_tract_freqs_tocat

   shell:
   '''
   module load !{params.mod_bcftools}
   #Calculate the median Sprime allele frequency for each tract:
   for chr in {!{autosome_list}};
      do
      if [[ "!{has_chr_prefix}" == "true" ]]; then
         chrnum=${chr#chr};
      else
         chrnum=${chr};
      fi;
      bcftools query -f '%CHROM:%POS[\t%GT]\n' -r ${chr} -H -S !{qpop}_samples.tsv !{qpop}_chr${chrnum}_nomissinggenos.vcf.gz | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeArchaicAF.awk -v "spop=!{pop}" -v "pop=!{qpop}" -v "all=1" !{pop}_chr${chrnum}_Sprime.score - | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractMedianAF.awk
   done | \
      gzip -9 > !{pop}_Sprime_!{qpop}_tract_freqs.tsv.gz
   '''
}

process sprime_genes {
   tag "${pop}"

   cpus params.sprimegenes_cpus
   memory { params.sprimegenes_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.sprimegenes_timeout }
   queue { task.attempt >= 2 ? params.sprimegenes_long_queue : params.sprimegenes_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

//   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.tcsv.gz'

   input:
   tuple val(pop), path("${pop}_Sprime_${pop}_tract_freqs.tsv.gz") from Sprime_tract_freqs.filter({ it[0] == it[1] }).map({ [it[0], it[2]] })
   path annotation_gff

   output:
   tuple val(pop), path("${pop}_Sprime_gene_lists.tcsv.gz"), path("${pop}_Sprime_gene_name_lists.tcsv.gz") into Sprime_gene_lists

   shell:
   '''
   module load !{params.mod_bedtools}
   #Identify genes overlapping Sprime tracts:
   bedtools intersect \
      -a <(gzip -dc !{pop}_Sprime_!{pop}_tract_freqs.tsv.gz | \
         awk 'BEGIN{FS="\t";OFS=FS;}{print $1, $2, $3, "TractID="$4";AAF="$6, ".", "+";}' | \
         sort -k1,1 -k2,2n -k3,3n) \
      -b <(gzip -dc !{annotation_gff} | \
         awk -v "keepchr=!{has_chr_prefix}" 'BEGIN{FS="\t";OFS=FS;}/^#/{print;}!/^#/{if (keepchr != "true") {sub("chr", "", $1);}; print $0;}') \
      -wao | \
      tee >(!{projectDir}/HumanPopGenScripts/Sprime/SprimeTractGeneList.awk -v "trim=1" -v "tag=gene_name" | \
         sort -k1,1V -k2,2n -k3,3n | \
         gzip -9 > !{pop}_Sprime_gene_name_lists.tcsv.gz) | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractGeneList.awk -v "trim=1" -v "tag=ID" | \
      sort -k1,1V -k2,2n -k3,3n | \
      gzip -9 > !{pop}_Sprime_gene_lists.tcsv.gz
   '''
}

process cat_outs {
   cpus params.catouts_cpus
   memory { params.catouts_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.catouts_timeout }
   queue { task.attempt >= 2 ? params.catouts_long_queue : params.catouts_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.t{sv,csv}.gz'

   input:
   path matchrates from Sprime_matchrates.collectFile() { [ "matchrate_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path projections from Sprime_project_BEDs.collectFile() { [ "projection_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path tractfreqs from Sprime_tract_freqs_tocat.map({ it[2] }).collectFile() { [ "tractfreq_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path genelists from Sprime_gene_lists.map({ it[2] }).collectFile() { [ "genelist_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path targetpops from sprime_target_pops.collectFile(name: 'Sprime_target_populations.txt', newLine: true, sort: true)
   path metadata
   path excluded_samples

   output:
   tuple path("${params.run_name}_perPop_Sprime_autosomal_match_rates.tsv.gz"), path("${params.run_name}_Sprime_perChrom_perIndiv_perPop_tract_lengths.tsv.gz"), path("${params.run_name}_Sprime_allPopPairs_tract_freqs.tsv.gz"), path("${params.run_name}_Sprime_gene_name_lists.tcsv.gz") into catouts_outputs

   shell:
   '''
   #Symlink the many input files to avoid SLURM SBATCH script size limits:
   cat matchrate_paths.tsv projection_paths.tsv tractfreq_paths.tsv genelist_paths.tsv | \
      while IFS=$'\t' read -a a;
         do
         ln -s ${a[2]} ${a[1]};
      done
   #Concatenate the match rate files across populations:
   header=1
   while IFS=$'\t' read p;
      do
      gzip -dc ${p}_autosomes_Sprime_match_rates.tsv.gz | \
         awk -v "header=${header}" 'BEGIN{FS="\t";OFS=FS;}NR==1&&header{print;}NR>1{print;}'
      header=0;
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_perPop_Sprime_autosomal_match_rates.tsv.gz
   unset header
   #Calculate total per-individual tract length and
   # concatenate files across populations:
   header=1
   while IFS=$'\t' read p;
      do
      !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=${p}" !{metadata} | \
         !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBEDtoLengths.awk -v "header=${header}" -v "pop=${p}" - ${p}_Sprime_tracts_perSample.bed
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_perChrom_perIndiv_perPop_tract_lengths.tsv.gz
   unset header
   #Concatenate tract frequency files and add a header:
   header=1
   while IFS=$'\t' read p;
      do
      if [[ "${header}" == "1" ]]; then
         printf "Chromosome\tStart\tEnd\tTractID\tQueryPop\tMedianAF\tMinAF\tMaxAF\tNumSites\n"
         header=""
      fi
      while IFS=$'\t' read q;
         do
         gzip -dc ${p}_Sprime_${q}_tract_freqs.tsv.gz | \
            sort -k1,1V -k2,2n -k3,3n;
      done < !{targetpops};
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_allPopPairs_tract_freqs.tsv.gz
   #Concatenate the gene name lists across populations:
   header=1
   while IFS=$'\t' read p;
      do
      if [[ "${header}" == "1" ]]; then
         printf "Chromosome\tStart\tEnd\tTractID\tOverlappingGenes\n"
         header=""
      fi
      gzip -dc ${p}_Sprime_gene_name_lists.tcsv.gz
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_gene_name_lists.tcsv.gz
   unset header
   '''
}

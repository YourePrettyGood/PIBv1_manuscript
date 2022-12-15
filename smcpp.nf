#!/usr/bin/env nextflow
/* Pipeline to infer single-population demographic histories with SMC++     *
 * Core steps:                                                              *
 *  Extraction of target population from main VCF for SMC++ ->              *
 *  Composing per-population masks from mosdepth's quantized BEDs |->       *
 *  smcpp vcf2smc with n distinguished lineages ->                          *
 *  smcpp estimate for each population and distinguished lineage ->         *
 *  smcpp plot --csv for each run of estimate                               *
 *  |-> bcftools +split to extract 4 samples individually from each pop ->  *
 *   generate_multihetsep.py to make MSMC2 inputs for each pop ->           *
 *   msmc2 on each pop                                                      *
 *  |-> CHIMP on each pop (from VCFs, no masks)                             */

//Default paths, globs, and regexes:
//Jointly genotyped VCFs:
params.vcf_glob = "${projectDir}/*.vcf.gz"
//Regex for extracting chromosome from VCF filename:
params.vcf_regex = ~/_chr(\p{Alnum}+)$/
//Glob for the per-individual BAMs and their indices:
params.bam_glob = "${projectDir}/BAMs/*_MD_IR_recal_filtered.ba{m,m.bai}"
//Regex for parsing the sample ID out of the BAM filename:
params.bam_regex = ~/^(.+)_MD_IR_recal_filtered$/
//Glob for the per-individual autosomal depth distributions from mosdepth
// generated by filter_VCF.nf:
params.depth_dist_glob = "${projectDir}/depthdists/*.mosdepth.global.dist.txt"
//Regex for parsing the sample ID out of the depth dist filename:
params.depth_regex = ~/^(.+)$/
//Sample ID column name in metadata file:
params.id_colname = "SampleID"
//Target group column name:
params.smcpp_target_colname = "AnalysisGroup"

//Regexes for CHIMP ref and anc FASTAs:
params.ref_bychrom_regex = ~/^(.+)$/
params.anc_bychrom_regex = ~/^.+_(.+?)$/

//Whether or not to run MSMC2:
params.run_msmc2 = false
//Whether or not to run CHIMP:
params.run_chimp = false

//Reference-related parameters for the pipeline:
params.ref_prefix = "/gpfs/gibbs/pi/tucci/pfr8/refs"
params.ref = "${params.ref_prefix}/1kGP/hs37d5/hs37d5.fa"
params.autosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"

//Set up the channels of per-chromosome jointly-genotyped and filtered VCFs and their indices:
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

//Set up the channel of populations for SMC++:
Channel
   .fromPath(params.smcpp_pops_file, checkIfExists: true)
   .ifEmpty { error "Unable to find populations file for SMC++: ${params.smcpp_pop_file}" }
   .splitCsv(sep: "\t")
   .map { it[0] }
   .tap { smcpp_pops }
   .tap { smcpp_pops_forfiltering }
   .subscribe { println "Added ${it} to smcpp_pops channel" }

//Set up a channel for the sample-to-population map to use for mask construction:
Channel
   .fromPath(params.metadata_file, checkIfExists: true)
   .ifEmpty { error "Unable to find metadata file: ${params.metadata_file}" }
   .splitCsv(sep: "\t", header: true)
   .map { [ it[params.smcpp_target_colname], it[params.id_colname] ] }
   .combine(smcpp_pops_forfiltering, by: 0)
   .tap { samples_forbamfiltering }
   .tap { samples_fordistfiltering }
   .map { [ it[1], it[0] ] }
   .tap { pop_map }
   .tap { pop_map_msmc2 }
   .subscribe { "Added ${it[0]}->${it[1]} to pop_map channel" }

//Set up the channel for all of the per-sample BAMs and their indices:
Channel
   .fromFilePairs(params.bam_glob, checkIfExists: true) { file -> (file.getSimpleName() =~ params.bam_regex)[0][1] }
   .ifEmpty { error "Unable to find BAMs matching glob: ${params.bam_glob}" }
   .combine(samples_forbamfiltering.map({ it[1] }), by: 0)
   .tap { bams }
   .subscribe { println "Added ${it[0]} to bams channel" }
//And a channel for the mosdepth *.global.dist.txt files from filter_VCF.nf:
Channel
   .fromPath(params.depth_dist_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find depth distribution files matching glob: ${params.depth_dist_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.depth_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract sample ID from depth distribution filename: ${params.depth_regex}" }
   .combine(samples_fordistfiltering.map({ it[1] }), by: 0)
   .tap { depth_dists }
   .subscribe { println "Added ${it[0]} (${it[1]}) to depth_dists channel" }

num_autosomes = params.autosomes.tokenize(',').size()

//Set up a value channel of the autosomes:
Channel
   .fromList(params.autosomes.tokenize(','))
   .tap { autosomes }
   .tap { autosomes_msmc2 }

//Set up the file channels for the metadata file:
metadata = file(params.metadata_file, checkIfExists: true)

//Set up the file channels for the ref and its various index components:
//Inspired by the IARC alignment-nf pipeline
//amb, ann, bwt, pac, and sa are all generated by bwa index, and alt is used by bwa when available for defining ALT mappings
//fai is generated by samtools faidx, and dict is generated by Picard and used by GATK
ref = file(params.ref, checkIfExists: true)
ref_dict = file(params.ref.replaceFirst("[.]fn?a(sta)?([.]gz)?", ".dict"), checkIfExists: true)
ref_fai = file(params.ref+'.fai', checkIfExists: true)

//A couple extra channels if CHIMP is to be run:
Channel
   .fromPath(params.ref_bychrom_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find per-chromosome reference FASTAs for CHIMP matching glob: ${params.ref_bychrom_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.ref_bychrom_regex)[0][1], a] }
   .filter { params.autosomes.tokenize(",").contains(it[0]) }
   .tap { refs }
   .subscribe { println "Added ${it[0]} (${it[1]}) to refs channel for CHIMP" }

Channel
   .fromPath(params.anc_bychrom_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find per-chromosome ancestral state FASTAs for CHIMP matching glob: ${params.anc_bychrom_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.anc_bychrom_regex)[0][1], a] }
   .filter { params.autosomes.tokenize(",").contains(it[0]) }
   .tap { ancs }
   .subscribe { println "Added ${it[0]} (${it[1]}) to ancs channel for CHIMP" }

//Default parameter values:
//Filters to apply to the VCF:
params.includestr = ""
//MAPQ threshold for calculating depth with mosdepth for GATK CallableLoci-like thresholding on DP:
params.mapq_thresh = "30"
//Lower and upper bounds for depth when making the per-sample masks:
params.mindp = "10"
params.maxdp_quantile = "0.995"
//Minimum fraction of samples marked uncallable required to mask a site:
//0 results in taking the union of per-sample uncallable masks
//1 results in taking the intersect of per-sample uncallable masks
params.frac_uncallable = "0"
//Mutation rate to use for SMC++:
params.mutation_rate = '1.25e-8'
//Mutation rate to use for CHIMP:
params.mut_rate = '0.0000000125'
//Recombination rate to use for CHIMP:
params.rec_rate = '0.0000000125'
//Seed for PRNGs:
params.prng_seed = 42
//Number of distinguished lineages to choose for SMC++ composite likelihood:
params.num_dlineages = 10
//Calculate the number of input .smc.gz files per SMC++ estimate call:
//This is simply the number of autosomes times the number of distinguished
// lineages to use for the composite likelihood.
smcpp_chroms_times_dlineages = num_autosomes.multiply(params.num_dlineages)
//Number of samples to use per population for MSMC2:
params.msmc2_num_samples = 2
params.msmc2_num_haplotypes = params.msmc2_num_samples.multiply(2)
//Time pattern string for MSMC2:
//Default is 1*2+25*1+1*2+1*3
params.msmc2_time_pattern = '1*2+25*1+1*2+1*3'
//Sample sizes to use for CHIMP (comma-separated list of integers):
params.chimp_n_s = '2,5,10'

//Defaults for cpus, memory, and time for each process:
//Constructing per-population mask BEDs
params.mosdepth_cpus = 1
params.mosdepth_mem = 2
params.mosdepth_timeout = '12h'
params.popmask_cpus = 1
params.popmask_mem = 1
params.popmask_timeout = '12h'
//VCF subsetting for SMC++
//Memory in MB
params.smcppsubset_cpus = 1
params.smcppsubset_mem = 512
params.smcppsubset_timeout = '1h'
//smcpp vcf2smc
//Memory in MB
params.vcftosmc_cpus = 1
params.vcftosmc_mem = 1024
params.vcftosmc_timeout = '12h'
//smcpp estimate
params.smcpp_cpus = 8
params.smcpp_mem = 32
params.smcpp_timeout = '24h'
//smcpp plot
params.smcppplot_cpus = 1
params.smcppplot_mem = 1
params.smcppplot_timeout = '6h'
//generate_multihetsep.py
params.msmcprep_cpus = 1
params.msmcprep_mem = 4
params.msmcprep_timeout = '24h'
//msmc2
params.msmc2_cpus = params.msmc2_num_haplotypes.multiply(params.msmc2_num_haplotypes.minus(1)).intdiv(2)
params.msmc2_mem = 40
params.msmc2_timeout = '24h'
//CHIMP
params.chimp_cpus = 1
params.chimp_mem = 8
params.chimp_timeout = '48h'

//Preprocess the per-chromosome VCF channel to include the indices:
perchrom_vcfs
   .join(perchrom_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
   .filter({ params.autosomes.tokenize(',').contains(it[0]) })
   .set { perchrom_vcfs_subset }

process smcpp_vcf_subset {
   tag "${pop} chr${chrom}"

   cpus params.smcppsubset_cpus
   memory { params.smcppsubset_mem.plus(task.attempt.minus(1).multiply(512))+' MB' }
   time { task.attempt >= 2 ? '24h' : params.smcppsubset_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), val(chrom), path(input_vcf), path(input_tbi) from smcpp_pops.combine(perchrom_vcfs_subset)
   path metadata

   output:
   tuple path("bcftools_view_smcpp_${pop}_chr${chrom}.stderr"), path("bcftools_view_smcpp_${pop}_chr${chrom}.stdout") into smcpp_vcf_subset_logs
   tuple val(pop), val(chrom), path("${pop}_chr${chrom}.vcf.gz"), path("${pop}_chr${chrom}.vcf.gz.tbi") into smcpp_vcfs,msmc2_vcfs,chimp_vcfs

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.smcpp_target_colname}" -v "select=!{pop}" !{metadata} > !{pop}_samples.tsv
   bcftools view -S !{pop}_samples.tsv -Oz -o !{pop}_chr!{chrom}.vcf.gz !{input_vcf} 2> bcftools_view_smcpp_!{pop}_chr!{chrom}.stderr > bcftools_view_smcpp_!{pop}_chr!{chrom}.stdout
   tabix -f !{pop}_chr!{chrom}.vcf.gz
   '''
}

process mosdepth {
   tag "${sampleid}"

   cpus params.mosdepth_cpus
   memory { params.mosdepth_mem.plus(task.attempt.minus(1).multiply(4))+' GB' }
   time { task.attempt == 2 ? '72h' : params.mosdepth_timeout }
   errorStrategy { task.exitStatus in ([1]+(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/callable", mode: 'copy', pattern: '*.bed.gz'

   input:
   tuple val(sampleid), path(depthdist), path(bambai) from depth_dists.join(bams, by: 0, failOnMismatch: true, failOnDuplicate: true)
   path ref
   path ref_fai

   output:
   tuple path("mosdepth_quantize_${sampleid}.stderr"), path("mosdepth_quantize_${sampleid}.stdout") into mosdepth_logs
   tuple val(sampleid), path("${sampleid}.mask.bed.gz") into mosdepth_masks
   tuple val(sampleid), path("${sampleid}.callable.bed.gz") into mosdepth_callable

   shell:
   '''
   module load !{params.mod_mosdepth}
   maxdp=$(awk 'BEGIN{FS="\t";OFS=FS;maxdpq=1-!{params.maxdp_quantile};}$1=="total"&&$3<=maxdpq{dpthresh=$2;}END{print dpthresh;}' !{depthdist})
   MOSDEPTH_Q0=NO_DEPTH MOSDEPTH_Q1=LOW_DEPTH MOSDEPTH_Q2=CALLABLE MOSDEPTH_Q3=HIGH_DEPTH mosdepth -n --fast-mode --mapq !{params.mapq_thresh} --quantize 0:1:!{params.mindp}:${maxdp}: !{sampleid}_quantized !{bambai[0]} 2> mosdepth_quantize_!{sampleid}.stderr > mosdepth_quantize_!{sampleid}.stdout
   gzip -dc !{sampleid}_quantized.quantized.bed.gz | \
      awk 'BEGIN{FS="\t";OFS=FS;}$4!="CALLABLE"{print $1, $2, $3;}' | \
      gzip -9 > !{sampleid}.mask.bed.gz
   gzip -dc !{sampleid}_quantized.quantized.bed.gz | \
      awk 'BEGIN{FS="\t";OFS=FS;}$4=="CALLABLE"{print $1, $2, $3;}' | \
      gzip -9 > !{sampleid}.callable.bed.gz
   '''
}

process pop_mask {
   tag "${pop} chr${chr}"

   cpus params.popmask_cpus
   memory { params.popmask_mem.plus(task.attempt.minus(1).multiply(4))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.popmask_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/uncallable", mode: 'copy', pattern: '*.bed.gz'

   input:
   tuple val(pop), path(perindiv_masks), val(chr) from pop_map
                                                          .join(mosdepth_masks, by: 0, failOnMismatch: true, failOnDuplicate: true)
                                                          .map { [ it[1], it[2] ] }
                                                          .groupTuple(by: 0)
                                                          .combine(autosomes)

   output:
   tuple val(pop), val(chr), path("${pop}_chr${chr}_mask.bed.gz"), path("${pop}_chr${chr}_mask.bed.gz.tbi") into pop_masks

   shell:
   '''
   module load !{params.mod_htslib}
   module load !{params.mod_bedtools}
   n_masks=$(ls !{perindiv_masks} | wc -l)
   bedtools multiinter -i !{perindiv_masks} | \
      awk -v "chrom=!{chr}" -v "n=${n_masks}" -v "thresh=!{params.frac_uncallable}" 'BEGIN{FS="\t";OFS=FS;}$1==chrom{if ($4/n>=thresh) {print $1, $2, $3;};}' | \
      bgzip -c > !{pop}_chr!{chr}_mask.bed.gz
   tabix -f !{pop}_chr!{chr}_mask.bed.gz
   '''
}

process vcf_to_smc {
   tag "${pop} chr${chrom}"

   cpus params.vcftosmc_cpus
   memory { params.vcftosmc_mem.plus(task.attempt.minus(1).multiply(1024))+' MB' }
   time { task.attempt >= 2 ? '24h' : params.vcftosmc_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), val(chrom), path("${pop}_chr${chrom}.vcf.gz"), path("${pop}_chr${chrom}.vcf.gz.tbi"), path(smcpp_mask_bed), path(smcpp_mask_bed_tbi) from smcpp_vcfs.join(pop_masks, by: [0,1], failOnMismatch: true, failOnDuplicate: true)

   output:
   tuple path("smcpp_vcf2smc_${pop}_chr${chrom}_d*.stderr"), path("smcpp_vcf2smc_${pop}_chr${chrom}_d*.stdout") into vcf_to_smc_logs
   tuple val(pop), path("${pop}_chr${chrom}_d*.smc.gz") into smcpp_inputs

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_smcpp}
   #Function for seeding PRNGs in utilities like shuf, shred, and sort:
   #See https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html
   seed_coreutils_PRNG()
   {
      seed="$1"
      openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt < /dev/zero 2> /dev/null
   }
   #Determine the number of samples in the input VCF:
   bcftools query -l !{pop}_chr!{chrom}.vcf.gz > !{pop}_samples.tsv
   sampleid_max=$(awk 'END{print NR-1;}' !{pop}_samples.tsv)
   shuf --random-source=<(seed_coreutils_PRNG !{params.prng_seed}) -i 0-${sampleid_max} -n!{params.num_dlineages} | \
   while read aid;
      do
      #Since we're working with unphased data, the distinguished lineages must come
      # from the same individual, so bid==aid:
      bid="${aid}"
      #Compose the distinguished lineage string and the population sample list:
      dlineage=""
      samplelist="!{pop}:"
      i=0
      while read sampleid <&3;
         do
         if [[ "${i}" -gt "0" ]]; then
            samplelist="${samplelist},"
         fi
         samplelist="${samplelist}${sampleid}"
         if [[ "${i}" -eq "${aid}" ]]; then
            dlineage="${sampleid}"
         fi
         if [[ "${i}" -eq "${bid}" ]]; then
            dlineage="${dlineage} ${sampleid}"
         fi
         i=$((i+1))
      done 3< !{pop}_samples.tsv
      smc++ vcf2smc -m !{smcpp_mask_bed} -d ${dlineage} !{pop}_chr!{chrom}.vcf.gz !{pop}_chr!{chrom}_d${aid}_${bid}.smc.gz !{chrom} ${samplelist} 2> smcpp_vcf2smc_!{pop}_chr!{chrom}_d${aid}_${bid}.stderr > smcpp_vcf2smc_!{pop}_chr!{chrom}_d${aid}_${bid}.stdout
   done
   '''
}

process smcpp_estimate {
   tag "${pop}"

   cpus params.smcpp_cpus
   memory { params.smcpp_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.smcpp_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '.debug.txt', saveAs: { 'smcpp_'+pop+it }
   publishDir path: "${params.output_dir}/smcpp", mode: 'copy', pattern: 'model.final.json', saveAs: { pop+'.'+it }
   publishDir path: "${params.output_dir}/smcpp", mode: 'copy', pattern: '.model.*.json', saveAs: { pop+it }

   input:
   tuple val(pop), path(smcinputpaths) from smcpp_inputs.flatMap({ it[1] })
                                                        .map({ [ (it.getSimpleName() =~ ~/^(.+?)_/)[0][1], it ] })
                                                        .collectFile() { [ "${it[0]}.txt", it[1].getSimpleName()+'\t'+it[1].getName()+'\t'+it[1]+'\n' ] }
                                                        .map({ [ it.getSimpleName(), it ] })
//   tuple val(pop), path(smcinputs) from smcpp_inputs.groupTuple(by: 0, size: smcpp_chroms_times_dlineages)

   output:
   tuple path("smcpp_estimate_${pop}.stderr"), path("smcpp_estimate_${pop}.stdout"), path(".debug.txt") into smcpp_estimate_logs
   tuple val(pop), path("model.final.json") into smcpp_final_models
   tuple val(pop), path(".model.iter*.json") into smcpp_intermediate_models

   shell:
   '''
   module load !{params.mod_smcpp}
   #Symlink the many input files from the current population !{pop}:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[2]} ${a[1]};
   done < !{smcinputpaths}
   #Now run SMC++ estimate on all of the inputs together:
   smc++ estimate --seed !{params.prng_seed} --cores !{task.cpus} -o . !{params.mutation_rate} *.smc.gz 2> smcpp_estimate_!{pop}.stderr > smcpp_estimate_!{pop}.stdout
   '''
}

process smcpp_plot {
   tag "${pop}"

   cpus params.smcppplot_cpus
   memory { params.smcppplot_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.smcppplot_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/smcpp", mode: 'copy', pattern: '*.{csv,png}'

   input:
   tuple val(pop), path("model.final.json"), path("*") from smcpp_final_models.join(smcpp_intermediate_models, by: 0, failOnDuplicate: true, failOnMismatch: true)

   output:
   tuple path("smcpp_plot_${pop}.stderr"), path("smcpp_plot_${pop}.stdout") into smcpp_plot_logs
   tuple val(pop), path("${pop}_all_models.png") into smcpp_perpop_plots
   tuple val(pop), path("${pop}_all_models.csv") into smcpp_perpop_csvs

   shell:
   '''
   module load !{params.mod_smcpp}
   smc++ plot --csv !{pop}_all_models.png model.final.json .model.iter*.json 2> smcpp_plot_!{pop}.stderr > smcpp_plot_!{pop}.stdout
   '''
}

process msmc2_prep {
   tag "${pop} chr${chrom}"

   cpus params.msmcprep_cpus
   memory { params.msmcprep_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.msmcprep_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   when: params.run_msmc2

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), val(chrom), path("*"), path("*"), path("*") from pop_map_msmc2
      .join(mosdepth_callable, by: 0, failOnMismatch: true, failOnDuplicate: true)
      .map { [ it[1], it[2] ] }
      .groupTuple(by: 0)
      .combine(autosomes_msmc2)
      .map { [ it[0], it[2], it[1] ] }
      .join(msmc2_vcfs, by: [0, 1], failOnMismatch: true, failOnDuplicate: true)

   output:
   tuple path("bcftools_view_bSNPs_msmc2_${pop}_samples_chr${chrom}.stderr"), path("bcftools_split_msmc2_${pop}_samples_chr${chrom}.stderr"), path("bcftools_split_msmc2_${pop}_samples_chr${chrom}.stdout"), path("generate_multihetsep_${pop}_chr${chrom}.stderr") into msmc2_prep_logs
   tuple val(pop), path("${pop}_chr${chrom}_multihetsep.txt") into msmc2_inputs

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_msmctools}
   #Function for seeding PRNGs in utilities like shuf, shred, and sort:
   #See https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html
   seed_coreutils_PRNG()
   {
      seed="$1"
      openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt < /dev/zero 2> /dev/null
   }
   #Get the list of samples to retain:
   bcftools query -l !{pop}_chr!{chrom}.vcf.gz | \
      shuf --random-source=<(seed_coreutils_PRNG !{params.prng_seed}) -n!{params.msmc2_num_samples} > samples_to_use.txt
   #Split the VCF into per-individual VCFs of these samples:
   bcftools view -Ou -v snps -m 2 -M 2 !{pop}_chr!{chrom}.vcf.gz 2> bcftools_view_bSNPs_msmc2_!{pop}_samples_chr!{chrom}.stderr | \
      bcftools +split -S samples_to_use.txt -Oz -o . - 2> bcftools_split_msmc2_!{pop}_samples_chr!{chrom}.stderr > bcftools_split_msmc2_!{pop}_samples_chr!{chrom}.stdout
   #Now generate the multihetsep files using the per-sample callable BEDs and VCFs:
   maskargs=""
   vcfargs=""
   while read sampleid;
      do
      maskargs="${maskargs} --mask=${sampleid}.callable.bed.gz"
      vcfargs="${vcfargs} ${sampleid}.vcf.gz"
   done < samples_to_use.txt
   generate_multihetsep.py ${maskargs} ${vcfargs} > !{pop}_chr!{chrom}_multihetsep.txt 2> generate_multihetsep_!{pop}_chr!{chrom}.stderr
   '''
}

process msmc2 {
   tag "${pop}"

   cpus params.msmc2_cpus
   memory { params.msmc2_mem.plus(task.attempt.minus(1).multiply(32))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.msmc2_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   when: params.run_msmc2

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.log'
   publishDir path: "${params.output_dir}/msmc2", mode: 'copy', pattern: '*.msmc2.*.txt'

   input:
   tuple val(pop), path(multihetseps) from msmc2_inputs.groupTuple(by: 0)

   output:
   tuple path("msmc2_${pop}.stderr"), path("msmc2_${pop}.stdout") into msmc2_logs
   tuple val(pop), path("${params.run_name}_${pop}_n${params.msmc2_num_haplotypes}.msmc2.loop.txt"), path("${params.run_name}_${pop}_n${params.msmc2_num_haplotypes}.msmc2.final.txt") into msmc2_output

   shell:
   '''
   module load !{params.mod_msmc2}
   haplist=$(seq -s, 0 $((!{params.msmc2_num_haplotypes}-1)))
   msmc2 -t !{task.cpus} -p "!{params.msmc2_time_pattern}" -o !{params.run_name}_!{pop}_n!{params.msmc2_num_haplotypes}.msmc2 -I ${haplist} !{multihetseps} 2> msmc2_!{pop}.stderr > msmc2_!{pop}.stdout
   '''
}

process chimp {
   tag "${pop}"

   cpus params.chimp_cpus
   memory { params.chimp_mem.plus(1).plus(task.attempt.minus(1).multiply(32))+' GB' }
   time { task.attempt >= 2 ? '120h' : params.chimp_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   when: params.run_chimp

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/CHIMP", mode: 'copy', pattern: '*.{csv,param}'

   input:
   tuple val(pop), path(vcfs), path(tbis), path(refs), path(ancs) from chimp_vcfs.map({ [ it[1], it[0], it[2], it[3] ] })
      .combine(refs.join(ancs, by: 0, failOnDuplicate: true, failOnMismatch: true), by: 0)
      .map({ [ it[1], it[2], it[3], it[4], it[5] ] })
      .groupTuple(by: 0)

   output:
   tuple path("CHIMP_autosomes_${pop}.stderr"), path("CHIMP_autosomes_${pop}.stdout") into chimp_logs
   tuple val(pop), path("CHIMP_autosomes_${pop}.csv"), path("CHIMP_autosomes_${pop}.param") into chimp_outputs

   shell:
   ref_list = refs.join("\n")
   anc_list = ancs.join("\n")
   vcf_list = vcfs.join("\n")
   chimp_retry_mem = params.chimp_mem.plus(task.attempt.minus(1).multiply(32))
   '''
   module load !{params.mod_chimp}
   #Generate sorted comma-separated lists for the per-chromosome inputs:
   ref_cslist=$(printf "!{ref_list}\n" | sort -k1,1V | head -c-1 | tr "\n" ",")
   anc_cslist=$(printf "!{anc_list}\n" | sort -k1,1V | head -c-1 | tr "\n" ",")
   vcf_cslist=$(printf "!{vcf_list}\n" | sort -k1,1V | head -c-1 | tr "\n" ",")
   #Run CHIMP with all the VCFs and corresponding ref and anc FASTAs:
   params="--rec_rate=!{params.rec_rate} --mut_rate=!{params.mut_rate} --base_n=!{params.chimp_n_s}"
   params="${params} --ref_list ${ref_cslist} --anc_list ${anc_cslist} --vcf_list ${vcf_cslist}"
   java -Xmx!{chimp_retry_mem}g -jar ${CHIMP} ${params} --out_file CHIMP_autosomes_!{pop} 2> CHIMP_autosomes_!{pop}.stderr > CHIMP_autosomes_!{pop}.stdout
   '''
}

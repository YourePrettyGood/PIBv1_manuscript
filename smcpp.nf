#!/usr/bin/env nextflow
/* Pipeline to infer single-population demographic histories with SMC++     *
 * Core steps:                                                              *
 *  Extraction of target population from main VCF for SMC++ ->              *
 *  Composing per-population masks from mosdepth's quantized BEDs ->        *
 *  smcpp vcf2smc with n distinguished lineages |->                         *
 *  smcpp estimate for each population and distinguished lineage ->         *
 *  smcpp plot --csv for each run of estimate                               *
 *  |> bootstrap vcf2smc outputs ->                                         *
 *   smcpp estimate each bootstrap ->                                       *
 *   smcpp plot each bootstrap                                              *
 * => rename and concatenate all results                                    */

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

//Detection mechanism for "chr" prefixes of autosomes:
has_chr_prefix = params.autosomes.count('chr') > 0
autosome_list = params.autosomes
autosome_num_list = has_chr_prefix ? params.autosomes.replaceAll('chr', '') : autosome_list
Map autosome_map = [params.autosomes.replaceAll('chr', '').tokenize(','), params.autosomes.tokenize(',')]
   .transpose()
   .collectEntries()

num_autosomes = params.autosomes.tokenize(',').size()

//Set up a value channel of the autosomes without chr prefix if detected:
Channel
   .fromList(autosome_num_list.tokenize(','))
   .tap { autosomes }

//Set up an optional file channel for the mappability mask for MSMC2:
//params.mappability_mask = "NOMASK"
//mappability_mask = file(params.mappability_mask)

//Set up the file channels for the metadata file:
metadata = file(params.metadata_file, checkIfExists: true)

//Set up the file channels for the ref and its various index components:
//Inspired by the IARC alignment-nf pipeline
//amb, ann, bwt, pac, and sa are all generated by bwa index, and alt is used by bwa when available for defining ALT mappings
//fai is generated by samtools faidx, and dict is generated by Picard and used by GATK
ref = file(params.ref, checkIfExists: true)
ref_dict = file(params.ref.replaceFirst("[.]fn?a(sta)?([.]gz)?", ".dict"), checkIfExists: true)
ref_fai = file(params.ref+'.fai', checkIfExists: true)

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
//Seed for PRNGs:
params.prng_seed = 42
//Number of distinguished lineages to choose for SMC++ composite likelihood:
params.num_dlineages = 10
//Calculate the number of input .smc.gz files per SMC++ estimate call:
//This is simply the number of autosomes times the number of distinguished
// lineages to use for the composite likelihood.
smcpp_chroms_times_dlineages = num_autosomes.multiply(params.num_dlineages)
/*//Number of bootstraps to perform:
params.num_bootstraps = 0
//Bootstrap block size (bp):
params.bootstrap_block_size = 5000000
//Number of blocks per bootstrapped chromosome:
//hs37d5 has ~2881 Mbp in the 22 autosomes, so int(2881/22)=130, and int(130/5)=26
params.num_blocks_per_chrom = 26
//Number of chromosomes per bootstrap:
params.num_chroms_per_bootstrap = 22
//PRNG seed for bootstrapping:
params.bootstrap_prng_seed = 42

//Set up value channel for the bootstraps:
if (params.num_bootstraps > 0) {
   Channel
      .of(1..params.num_bootstraps)
      .tap { bootstraps }
      .subscribe { println "Adding bootstrap ${it} to bootstraps channel" }
} else {
   bootstraps = Channel.empty()
}*/

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
/*//bootstrap
params.bootstrap_cpus = 1
params.bootstrap_mem = 1
params.bootstrap_timeout = '2h'*/

//Preprocess the per-chromosome VCF channel to include the indices:
perchrom_vcfs
   .join(perchrom_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
   .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
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
   tuple val(pop), val(chrom), path("${pop}_chr${chrom}.vcf.gz"), path("${pop}_chr${chrom}.vcf.gz.tbi") into smcpp_vcfs

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
   tuple val(pop), path(perindiv_masks), val(chrom) from pop_map
                                                          .join(mosdepth_masks, by: 0, failOnMismatch: true, failOnDuplicate: true)
                                                          .map { [ it[1], it[2] ] }
                                                          .groupTuple(by: 0)
                                                          .combine(autosomes)

   output:
   tuple val(pop), val(chrom), path("${pop}_chr${chrom}_mask.bed.gz"), path("${pop}_chr${chrom}_mask.bed.gz.tbi") into pop_masks

   shell:
   chrid = autosome_map[chrom]
   '''
   module load !{params.mod_htslib}
   module load !{params.mod_bedtools}
   n_masks=$(ls !{perindiv_masks} | wc -l)
   bedtools multiinter -i !{perindiv_masks} | \
      awk -v "chrom=!{chrid}" -v "n=${n_masks}" -v "thresh=!{params.frac_uncallable}" 'BEGIN{FS="\t";OFS=FS;}$1==chrom{if ($4/n>=thresh) {print $1, $2, $3;};}' | \
      bgzip -c > !{pop}_chr!{chrom}_mask.bed.gz
   tabix -f !{pop}_chr!{chrom}_mask.bed.gz
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
   chrid = autosome_map[chrom]
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
   n_samples=$(awk 'END{print NR;}' !{pop}_samples.tsv)
   n_draws=!{params.num_dlineages}
   #Adjust the number of distinguished lineages considered if there
   # aren't enough samples:
   if [[ "${n_samples}" -lt "${n_draws}" ]]; then
      n_draws=${n_samples}
      echo "Warning: !{pop} only has ${n_samples} individuals, so cannot draw !{params.num_dlineages} distinguished lineages"
      echo "Drawing ${n_draws} distinguished lineages instead"
   fi
   sampleid_max=$(awk 'END{print NR-1;}' !{pop}_samples.tsv)
   shuf --random-source=<(seed_coreutils_PRNG !{params.prng_seed}) -i 0-${sampleid_max} -n${n_draws} | \\
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
      smc++ vcf2smc -m !{smcpp_mask_bed} -d ${dlineage} !{pop}_chr!{chrom}.vcf.gz !{pop}_chr!{chrom}_d${aid}_${bid}.smc.gz !{chrid} ${samplelist} 2> smcpp_vcf2smc_!{pop}_chr!{chrom}_d${aid}_${bid}.stderr > smcpp_vcf2smc_!{pop}_chr!{chrom}_d${aid}_${bid}.stdout
   done
   '''
}

/*if (params.num_bootstraps > 0) {
   smcpp_inputs
      .tap { smcpp_inputs_for_bootstrap }
} else {
   smcpp_inputs_for_bootstrap = Channel.empty()
}*/

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

/*process bootstrap {
   tag "${pop}"

   cpus params.bootstrap_cpus
   memory { params.bootstrap_mem.plus(task.attempt.minus(1).multiply(4))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.bootstrap_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(pop), path(smcinputpaths), val(bootstrap) from smcpp_inputs_for_bootstraps.flatMap({ it[1] })
                                                               .map({ [ (it.getSimpleName() =~ ~/^(.+?)_/)[0][1], it ] })
                                                               .collectFile() { [ "${it[0]}.txt", it[1].getSimpleName()+'\t'+it[1].getName()+'\t'+it[1]+'\n' ] }
                                                               .map({ [ it.getSimpleName(), it ] })
                                                               .combine(bootstraps)

   output:
   tuple val(pop), val(bootstrap), path("${pop}_b${bootstrap}_chr${chrom}_d*.smc.gz") into bootstrap_smcpp_inputs

   shell:
   '''
   #Symlink the many input files from the current population !{pop}:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[2]} ${a[1]};
   done < !{smcinputpaths}
   //Bootstrap separately for each distinguished lineage:
   awk 'BEGIN{FS="\t";OFS=FS;}{n=split(\$1, a, "_");print a[n];}' !{smcinputpaths} | \\
      while read d;
         do
         //Determine bootstrap chunk boundaries:
         gzip -dc *_${d}.smc.gz | \\
            !{projectDir}/HumanPopGenScripts/SMCpp/defineBootstrapBlocks.awk -v "blocksize=!{params.bootstrap_block_size}" > bootstrap_block_boundaries_${d}.txt
         //Perform the bootstrap:
         gzip -dc *_${d}.smc.gz | \\
            !{projectDir}/HumanPopGenScripts/SMCpp/bootstrapBlocks.awk -v "blocksperchrom=!{params.num_blocks_per_chrom}" -v "nchroms=!{params.num_chroms_per_bootstrap}" -v "prefix=!{pop}_b" -v "suffix=_${d}.smc.gz" -v "seed=!{params.bootstrap_prng_seed}" bootstrap_block_boundaries_${d}.txt -
      done
   '''
}

process bootstrap_estimate {
   tag "${pop}"

   cpus params.smcpp_cpus
   memory { params.smcpp_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.smcpp_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '.debug.txt', saveAs: { 'bootstrap_'+pop+it }
   publishDir path: "${params.output_dir}/bootstraps", mode: 'copy', pattern: 'model.final.json', saveAs: { pop+'.'+it }
   publishDir path: "${params.output_dir}/bootstraps", mode: 'copy', pattern: '.model.*.json', saveAs: { pop+it }

   input:
   tuple val(pop), path(smcinputpaths) from bootstrap_smcpp_inputs.flatMap({ it[2] })
                                               .map({ [ (it.getSimpleName() =~ ~/^(.+?_b[0-9]+)_/)[0][1], it ] })
                                               .collectFile() { [ "${it[0]}.txt", it[1].getSimpleName()+'\t'+it[1].getName()+'\t'+it[1]+'\n' ] }
                                               .map({ [ it.getSimpleName(), it ] })

   output:
   tuple path("bootstrap_estimate_${pop}.stderr"), path("bootstrap_estimate_${pop}.stdout"), path(".debug.txt") into bootstrap_estimate_logs
   tuple val(pop), path("model.final.json") into bootstrap_final_models
   tuple val(pop), path(".model.iter*.json") into bootstrap_intermediate_models

   shell:
   '''
   module load !{params.mod_smcpp}
   #Symlink the many input files from the current population !{pop}:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[2]} ${a[1]};
   done < !{smcinputpaths}
   #Now run SMC++ estimate on all of the inputs together:
   smc++ estimate --seed !{params.prng_seed} --cores !{task.cpus} -o . !{params.mutation_rate} *.smc.gz 2> bootstrap_estimate_!{pop}.stderr > bootstrap_estimate_!{pop}.stdout
   '''
}

process boot_plot {
   tag "${pop}"

   cpus params.smcppplot_cpus
   memory { params.smcppplot_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '48h' : params.smcppplot_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/bootstraps", mode: 'copy', pattern: '*.{csv,png}'

   input:
   tuple val(pop), path("model.final.json"), path("*") from bootstrap_final_models.join(bootstrap_intermediate_models, by: 0, failOnDuplicate: true, failOnMismatch: true)

   output:
   tuple path("bootstrap_plot_${pop}.stderr"), path("bootstrap_plot_${pop}.stdout") into bootstrap_plot_logs
   tuple val(pop), path("${pop}_all_models.png") into bootstrap_perpop_plots
   tuple val(pop), path("${pop}_all_models.csv") into bootstrap_perpop_csvs

   shell:
   '''
   module load !{params.mod_smcpp}
   smc++ plot --csv !{pop}_all_models.png model.final.json .model.iter*.json 2> bootstrap_plot_!{pop}.stderr > bootstrap_plot_!{pop}.stdout
   '''
}*/

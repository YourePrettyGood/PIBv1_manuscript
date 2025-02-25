#!/usr/bin/env nextflow
/* Pipeline for running hmmix for archaic introgression   *
 * Core steps:                                            *
 *  vcftobcf -> create_outgroup -> cat_ogs ->             *
 *  mutation_rate -> create_ingroup -> cat_igs ->         *
 *  train -> decode                                       *
 * Options include:                                       *
 *  haploid training and decoding (params.hmmix_haploid)  *
 *  annotate tracts with archaics (params.annotate_arc)   */

nextflow.enable.dsl = 2

//Default paths, globs, and regexes:
//Modern sample VCF glob:
params.modern_vcf_glob = "${projectDir}/modern_VCFs/*.vcf.gz"
//Regex to extract chromosome number from modern sample VCF:
params.modern_vcf_regex = ~/_chr(\p{Alnum}+)/
//Archaic sample VCF glob:
params.arc_vcf_glob = "${projectDir}/archaic_VCFs/*.vcf.gz"
//Regex to extract chromosome number from archaic sample VCF:
params.arc_vcf_regex = ~/_chr(\p{Alnum}+)/
//Ancestral allele FASTA glob:
params.anc_glob = "${projectDir}/AA_FASTAs/*.fa"
//Regex to extract chromosome number from ancestral allele FASTAs:
params.anc_regex = ~/_(\p{Alnum}+)$/
//Ref allele FASTA glob:
params.ref_glob = "${projectDir}/ref_FASTAs/*.fa"
//Regex to extract chromosome number from ref allele FASTAs:
params.ref_regex = ~/^chr(\p{Alnum}+)$/
//Also remember to set:
//Path to VCF mask (in BED format, positive mask so only intervals to be retained):
//params.vcf_mask = ""
//Path to JSON of ingroup and outgroup sample IDs required by hmmix:
//params.ind_json = ""

//Default parameters for hmmix:
//Train HMM on phased data to predict haplotype tracts?:
params.hmmix_haploid = false
//Annotate tracts with archaic allele sharing?:
params.annotate_arc = false
//Window size for estimating mutation rate (default: 1000000):
params.mutrate_window_size = "1000000"

//Defaults for cpus, memory, time, and queue for each process:
//Convert input VCFs to BCFs (vcftobcf):
// (GNU parallel -j4 took avg 127.0s for 1kGP 2-22, so 2,667s total with 8GB allocated by SLURM)
// (Archaics 1-22 took avg 140.3s, so 3,086.6s total)
// (.: 1h is almost surely safe, 2h is extra conservative, and 4GB is probably plenty)
params.vcftobcf_cpus = 1
params.vcftobcf_mem = 4
params.vcftobcf_mem_increment = 8
params.vcftobcf_timeout = '2h'
params.vcftobcf_retry_timeout = '8h'
params.vcftobcf_queue = 'day'
params.vcftobcf_retry_queue = 'day'
//Identify sites segregating in outgroup per-chromosome (outgroup):
params.outgroup_cpus = 1
params.outgroup_mem = 2
params.outgroup_mem_increment = 4
params.outgroup_timeout = '1h'
params.outgroup_retry_timeout = '6h'
params.outgroup_queue = 'day'
params.outgroup_retry_queue = 'day'
//Combine outgroup segregating sites across chromosomes (cat_ogs):
params.cat_ogs_cpus = 1
params.cat_ogs_mem = 1
params.cat_ogs_mem_increment = 4
params.cat_ogs_timeout = '1h'
params.cat_ogs_retry_timeout = '6h'
params.cat_ogs_queue = 'day'
params.cat_ogs_retry_queue = 'day'
//Estimate local mutation rate by density of subs in outgroup (mutrate):
params.mutrate_cpus = 1
params.mutrate_mem = 1
params.mutrate_mem_increment = 4
params.mutrate_timeout = '1h'
params.mutrate_retry_timeout = '6h'
params.mutrate_queue = 'day'
params.mutrate_retry_queue = 'day'
//Identify segsites in ingroup not found in outgroup per-chromosome (ingroup):
params.ingroup_cpus = 1
params.ingroup_mem = 4
params.ingroup_mem_increment = 8
params.ingroup_timeout = '1h'
params.ingroup_retry_timeout = '6h'
params.ingroup_queue = 'day'
params.ingroup_retry_queue = 'day'
//Combine ingroup segregating sites across chromosomes (cat_igs):
params.cat_igs_cpus = 1
params.cat_igs_mem = 1
params.cat_igs_mem_increment = 4
params.cat_igs_timeout = '1h'
params.cat_igs_retry_timeout = '6h'
params.cat_igs_queue = 'day'
params.cat_igs_retry_queue = 'day'
//Train the HMM for each individual (train):
// (<4m and <4GB per indiv on 1kGP)
params.train_cpus = 1
params.train_mem = 4
params.train_mem_increment = 8
params.train_timeout = '1h'
params.train_retry_timeout = '6h'
params.train_queue = 'day'
params.train_retry_queue = 'day'
//Perform posterior decoding for each individual (decode):
//Annotating with archaics increases runtime but not memory
// (<1m -> 8-9m, but <4GB for both on 1kGP)
params.decode_cpus = 1
params.decode_mem = 4
params.decode_mem_increment = 8
params.decode_timeout = '1h'
params.decode_retry_timeout = '6h'
params.decode_queue = 'day'
params.decode_retry_queue = 'day'

process vcftobcf {
   tag "${chrom} ${src}"

   cpus params.vcftobcf_cpus
   memory { params.vcftobcf_mem.plus(task.attempt.minus(1).multiply(params.vcftobcf_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.vcftobcf_retry_timeout : params.vcftobcf_timeout }
   queue { task.attempt >= 2 ? params.vcftobcf_retry_queue : params.vcftobcf_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(src), val(chrom), path(vcf), path(tbi)

   output:
   tuple path("bcftools_view_tobcf_${params.run_name}_${src}_${chrom}.stderr"), path("bcftools_view_tobcf_${params.run_name}_${src}_${chrom}.stdout"), path("bcftools_index_bcf_${params.run_name}_${src}_${chrom}.stderr"), path("bcftools_index_bcf_${params.run_name}_${src}_${chrom}.stdout"), emit: logs
   tuple val(src), val(chrom), path("${params.run_name}_${src}_${chrom}.bcf"), path("${params.run_name}_${src}_${chrom}.bcf.csi"), emit: bcf

   shell:
   '''
   module load !{params.mod_bcftools}
   bcftools view -Ob -o !{params.run_name}_!{src}_!{chrom}.bcf !{vcf} 2> bcftools_view_tobcf_!{params.run_name}_!{src}_!{chrom}.stderr > bcftools_view_tobcf_!{params.run_name}_!{src}_!{chrom}.stdout
   bcftools index !{params.run_name}_!{src}_!{chrom}.bcf 2> bcftools_index_bcf_!{params.run_name}_!{src}_!{chrom}.stderr > bcftools_index_bcf_!{params.run_name}_!{src}_!{chrom}.stdout
   '''
}

process outgroup {
   tag "${chrom}"

   cpus params.outgroup_cpus
   memory { params.outgroup_mem.plus(task.attempt.minus(1).multiply(params.outgroup_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.outgroup_retry_timeout : params.outgroup_timeout }
   queue { task.attempt >= 2 ? params.outgroup_retry_queue : params.outgroup_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(bcf), path(csi), path(anc), path(ref)
   path(mask)
   path(indjson)

   output:
   tuple path("hmmix_create_outgroup_${params.run_name}_${chrom}.stderr"), path("hmmix_create_outgroup_${params.run_name}_${chrom}.stdout"), emit: logs
   path("outgroup_${chrom}.txt"), emit: outgroup

   shell:
   '''
   module load !{params.mod_miniconda}
   conda activate hmmix
   hmmix create_outgroup -ind=!{indjson} -vcf=!{bcf} -weights=!{mask} -out=outgroup_!{chrom}.txt -ancestral=!{anc} -refgenome=!{ref} 2> hmmix_create_outgroup_!{params.run_name}_!{chrom}.stderr > hmmix_create_outgroup_!{params.run_name}_!{chrom}.stdout
   '''
}

process cat_ogs {
   cpus params.cat_ogs_cpus
   memory { params.cat_ogs_mem.plus(task.attempt.minus(1).multiply(params.cat_ogs_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.cat_ogs_retry_timeout : params.cat_ogs_timeout }
   queue { task.attempt >= 2 ? params.cat_ogs_retry_queue : params.cat_ogs_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   path(outgroups)

   output:
   path("outgroup.txt"), emit: outgroup

   shell:
   '''
   awk 'BEGIN{FS="\\t";OFS=FS;}NR==1{print;}FNR>1{print;}' !{outgroups} > outgroup.txt
   '''
}

process mutrate {
   cpus params.mutrate_cpus
   memory { params.mutrate_mem.plus(task.attempt.minus(1).multiply(params.mutrate_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.mutrate_retry_timeout : params.mutrate_timeout }
   queue { task.attempt >= 2 ? params.mutrate_retry_queue : params.mutrate_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   path(outgroup)
   path(mask)

   output:
   tuple path("hmmix_mutation_rate_${params.run_name}.stderr"), path("hmmix_mutation_rate_${params.run_name}.stdout"), emit: logs
   path("${params.run_name}_mutation_rate.bed"), emit: mutrate

   shell:
   '''
   module load !{params.mod_miniconda}
   conda activate hmmix
   hmmix mutation_rate -outgroup=!{outgroup} -weights=!{mask} -window_size=!{params.mutrate_window_size} -out !{params.run_name}_mutation_rate.bed 2> hmmix_mutation_rate_!{params.run_name}.stderr > hmmix_mutation_rate_!{params.run_name}.stdout
   '''
}

process ingroup {
   tag "${chrom}"

   cpus params.ingroup_cpus
   memory { params.ingroup_mem.plus(task.attempt.minus(1).multiply(params.ingroup_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.ingroup_retry_timeout : params.ingroup_timeout }
   queue { task.attempt >= 2 ? params.ingroup_retry_queue : params.ingroup_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(bcf), path(csi), path(anc), path(outgroup)
   path(mask)
   path(indjson)

   output:
   tuple path("hmmix_create_ingroup_${params.run_name}_${chrom}.stderr"), path("hmmix_create_ingroup_${params.run_name}_${chrom}.stdout"), emit: logs
   path("*_${chrom}_obs.txt"), emit: obs

   shell:
   '''
   module load !{params.mod_miniconda}
   conda activate hmmix
   hmmix create_ingroup -ind=!{indjson} -vcf=!{bcf} -weights=!{mask} -out=obs_!{chrom} -outgroup=!{outgroup} -ancestral=!{anc} 2> hmmix_create_ingroup_!{params.run_name}_!{chrom}.stderr > hmmix_create_ingroup_!{params.run_name}_!{chrom}.stdout
   for i in obs_!{chrom}.*.txt;
      do
      prefix=$(basename ${i} .txt);
      sampleid=${prefix#obs_!{chrom}.};
      mv ${i} ${sampleid}_!{chrom}_obs.txt;
   done
   '''
}

process cat_igs {
   tag "${id}"

   cpus params.cat_igs_cpus
   memory { params.cat_igs_mem.plus(task.attempt.minus(1).multiply(params.cat_igs_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.cat_igs_retry_timeout : params.cat_igs_timeout }
   queue { task.attempt >= 2 ? params.cat_igs_retry_queue : params.cat_igs_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(id), path(ingroups)

   output:
   tuple val(id), path("${id}_obs.txt"), emit: ingroup

   shell:
   '''
   awk 'BEGIN{FS="\\t";OFS=FS;}NR==1{print;}FNR>1{print;}' !{ingroups} > !{id}_obs.txt
   '''
}

process train {
   tag "${id}"

   cpus params.train_cpus
   memory { params.train_mem.plus(task.attempt.minus(1).multiply(params.train_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.train_retry_timeout : params.train_timeout }
   queue { task.attempt >= 2 ? params.train_retry_queue : params.train_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(id), path(obs), path(mutrate)
   path(mask)

   output:
   tuple path("hmmix_train_${params.run_name}_${id}${params.hap_arg}.stderr"), path("hmmix_train_${params.run_name}_${id}${params.hap_arg}.stdout"), emit: logs
   tuple val(id), path(obs), path("${params.run_name}_trained${params.hap_arg}.${id}.json"), path(mutrate), emit: hmmparams

   shell:
   '''
   module load !{params.mod_miniconda}
   conda activate hmmix
   hmmix train -obs=!{obs} -weights=!{mask} -mutrates=!{mutrate} !{params.hap_arg} -out=!{params.run_name}_trained!{params.hap_arg}.!{id}.json 2> hmmix_train_!{params.run_name}_!{id}!{params.hap_arg}.stderr > hmmix_train_!{params.run_name}_!{id}!{params.hap_arg}.stdout
   '''
}

process decode {
   tag "${id}"

   cpus params.decode_cpus
   memory { params.decode_mem.plus(task.attempt.minus(1).multiply(params.decode_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.decode_retry_timeout : params.decode_timeout }
   queue { task.attempt >= 2 ? params.decode_retry_queue : params.decode_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/hmmix", mode: 'copy', pattern: '*decoded*.txt'

   input:
   tuple val(id), path(obs), path(hmmparams), path(mutrate), path(arcbcfs), path(arccsis)
   path(mask)

   output:
   tuple path("hmmix_decode_${params.run_name}_${id}${params.hap_arg}${params.admixpop_arg}.stderr"), path("hmmix_decode_${params.run_name}_${id}${params.hap_arg}${params.admixpop_arg}.stdout"), emit: logs
   tuple val(id), path("${params.run_name}_${id}.decoded${params.admixpop_arg}.*.txt"), emit: calls

   shell:
   arcbcf_list = params.annotate_arc ? arcbcfs.join("\\n") : ""
   '''
   module load !{params.mod_miniconda}
   conda activate hmmix
   #Generate sorted comma-separated list of archaic BCFs if requested:
   arc="!{params.admixpop_arg}"
   if [[ -n "!{params.admixpop_arg}" ]]; then
      arcbcf_cslist=$(printf "!{arcbcf_list}\\n" | sort -k1,1V | head -c-1 | tr "\\n" ",");
      arc="${arc}=${arcbcf_cslist}";
   fi
   hmmix decode -obs=!{obs} -weights=!{mask} -mutrates=!{mutrate} -param=!{hmmparams} !{params.hap_arg} ${arc} -out=!{params.run_name}_!{id}.decoded!{params.admixpop_arg} 2> hmmix_decode_!{params.run_name}_!{id}!{params.hap_arg}!{params.admixpop_arg}.stderr > hmmix_decode_!{params.run_name}_!{id}!{params.hap_arg}!{params.admixpop_arg}.stdout
   '''
}

workflow {
   //Detection mechanism for "chr" prefixes of autosomes:
   has_chr_prefix = params.autosomes.count('chr') > 0
   autosome_list = params.autosomes
   autosome_num_list = has_chr_prefix ? params.autosomes.replaceAll('chr', '') : autosome_list
   num_autosomes = autosome_num_list.tokenize(',').size()

   //Set up the arguments for the options:
   params.hap_arg = params.hmmix_haploid ? "-haploid" : ""
   params.admixpop_arg = params.annotate_arc ? "-admixpop" : ""

   //Set up the channel of modern VCFs and indices:
   Channel
      .fromPath(params.modern_vcf_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find VCFs of modern samples matching glob: ${params.modern_vcf_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.modern_vcf_regex)[0][1], a ] }
      .ifEmpty { error "Regex failed to extract chromosome ID from modern sample VCF: ${params.modern_vcf_regex}" }
      .tap { modern_vcfs_only }
      .subscribe { println "Added ${it[1]} to modern_vcfs_only channel" }
   Channel
      .fromPath(params.modern_vcf_glob+'.tbi', checkIfExists: true)
      .ifEmpty { error "Unable to find VCF indices of modern samples matching glob: ${params.modern_vcf_glob}.tbi" }
      .map { a -> [ (a.getSimpleName() =~ params.modern_vcf_regex)[0][1], a ] }
      .ifEmpty { error "Regex failed to extract chromosome ID from modern sample VCF index: ${params.modern_vcf_regex}" }
      .tap { modern_tbis_only }
      .subscribe { println "Added ${it[1]} to modern_tbis_only channel" }
   modern_vcfs = modern_vcfs_only
      .join(modern_tbis_only, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
      .map { [ "modern", it[0], it[1], it[2] ] }

   if (params.annotate_arc) {
      //Set up the channel of archaic VCFs and indices:
      Channel
         .fromPath(params.arc_vcf_glob, checkIfExists: true)
         .ifEmpty { error "Unable to find VCFs of archaic samples matching glob: ${params.arc_vcf_glob}" }
         .map { a -> [ (a.getSimpleName() =~ params.arc_vcf_regex)[0][1], a ] }
         .ifEmpty { error "Regex failed to extract chromosome ID from archaic sample VCF: ${params.arc_vcf_regex}" }
         .tap { arc_vcfs_only }
         .subscribe { println "Added ${it[1]} to arc_vcfs_only channel" }
      Channel
         .fromPath(params.arc_vcf_glob+'.tbi', checkIfExists: true)
         .ifEmpty { error "Unable to find VCF indices of archaic samples matching glob: ${params.arc_vcf_glob}.tbi" }
         .map { a -> [ (a.getSimpleName() =~ params.arc_vcf_regex)[0][1], a ] }
         .ifEmpty { error "Regex failed to extract chromosome ID from archaic sample VCF index: ${params.arc_vcf_regex}" }
         .tap { arc_tbis_only }
         .subscribe { println "Added ${it[1]} to arc_tbis_only channel" }
      arc_vcfs = arc_vcfs_only
         .join(arc_tbis_only, by: 0, failOnDuplicate: true, failOnMismatch: true)
         .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
         .map { [ "archaic", it[0], it[1], it[2] ] }
   } else {
      arc_vcfs = Channel.empty()
   }

   //Set up the channel of ancestral allele FASTAs:
   Channel
      .fromPath(params.anc_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find ancestral allele FASTAs matching glob: ${params.anc_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.anc_regex)[0][1], a ] }
      .ifEmpty { error "Regex failed to extract chromosome ID from ancestral allele FASTA: ${params.anc_glob} =~ ${params.anc_regex}" }
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
      .tap { ancs_for_outgroup }
      .tap { ancs_for_ingroup }
      .subscribe { println "Added chr${it[0]} ancestral allele FASTA (${it[1]}) to ancs_* channels" }

   //Set up the channel of ref allele FASTAs:
   Channel
      .fromPath(params.ref_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find ref allele FASTAs matching glob: ${params.ref_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.ref_regex)[0][1], a ] }
      .ifEmpty { error "Regex failed to extract chromosome ID from ref allele FASTA: ${params.ref_glob} =~ ${params.ref_regex}" }
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
      .tap { refs }
      .subscribe { println "Added chr${it[0]} ref allele FASTA (${it[1]}) to refs channel" }

   //Set up the value channels of the mask BED file and the individuals JSON file:
   mask = file(params.vcf_mask, checkIfExists: true)
   indjson = file(params.ind_json, checkIfExists: true)

   //Convert input VCFs to BCFs:
   vcfs = modern_vcfs.mix(arc_vcfs)
   vcftobcf(vcfs)
   bcfs = vcftobcf.out.bcf
      .branch { src,chr,bcf,csi ->
         modern: src == "modern"
         arc: src == "archaic"
      }
   bcfs_for_outgroup = bcfs.modern
      .map { src,chr,bcf,csi -> [ chr, bcf, csi ] }
      .tap { bcfs_for_ingroup }
   arc_bcfs = bcfs.arc
      .ifEmpty({ [ "archaic", "1", file('no_arc.bcf'), file('no_arc.bcf.csi') ] })
      .map { src,chr,bcf,csi -> [ src, bcf, csi ] }
      .groupTuple(by: 0)
      .map { src,bcfs,csis -> [ bcfs, csis ] }

   //Generate the list of sites with outgroup variation:
   outgroup_inputs = bcfs_for_outgroup
      .join(ancs_for_outgroup, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .join(refs, by: 0, failOnDuplicate: true, failOnMismatch: true)
   outgroup(outgroup_inputs, mask, indjson)
   cat_ogs(outgroup.out.outgroup.collect())
   outgroup_for_mutrate = cat_ogs.out.outgroup
      .tap { outgroup_for_ingroup }

   //Estimate the mutation rate based on outgroup substitution density:
   mutrate(outgroup_for_mutrate, mask)

   //Generate the observations for each individual and concatenate across chroms:
   ingroup_inputs = bcfs_for_ingroup
      .join(ancs_for_ingroup, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .combine(outgroup_for_ingroup)
   ingroup(ingroup_inputs, mask, indjson)
   cat_igs_inputs = ingroup.out.obs
      .flatten()
      .map({ [ (it.getName() =~ ~/^(\w+)_\p{Alnum}+_obs.txt$/)[0][1], it ] })
      .groupTuple(by: 0, size: num_autosomes)
   cat_igs(cat_igs_inputs)

   //Train the HMM on each individual/pair of haplotypes using forward-backward:
   train_inputs = cat_igs.out.ingroup
      .combine(mutrate.out.mutrate)
   train(train_inputs, mask)
   //Now do the posterior decoding for each individual/pair of haplotypes:
   decode_inputs = train.out.hmmparams
      .combine(arc_bcfs)
   decode(decode_inputs, mask)
}


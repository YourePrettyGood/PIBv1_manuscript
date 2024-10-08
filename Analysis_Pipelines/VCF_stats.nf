#!/usr/bin/env nextflow
/* Pipeline for basic VCF stats                                          *
 * Steps:                                                                *
 *  bcftools stats                                                       *
 *  bcftools +smpl-stats                                                 *
 *  bcftools +trio-stats                                                 */

//Filtering variables:
params.filterstr = ""
params.includestr = ""

//Suffix for the current run's output logs:
params.run_suffix = "nofilters"

//Default paths, globs, and regexes:
params.vcf_glob = "${projectDir}/byChrom/*.vcf.gz"
//Regex to pull out the chromosome:
params.vcfchr_regex = ~/_chr(\p{Alnum}+)$/
//For the trio stats:
params.trio_ped = "${projectDir}/trio/PIBv1_trios.ped"
//For the array data for gtcheck:
params.array_vcf = "${projectDir}/array/TucciCombinedDatasetHumanOrigins_noreffix_flipswap_headercontigsfixed.vcf.gz"
params.array_sites = "${projectDir}/array/TucciCombinedDatasetHumanOrigins_sites.tsv"

//Reference-related parameters for the pipeline:
params.ref_prefix = "/gpfs/gibbs/pi/tucci/pfr8/refs"
params.ref = "${params.ref_prefix}/1kGP/hs37d5/hs37d5.fa"

//Set up the channel of VCFs
Channel
   .fromFilePairs(params.vcf_glob, size: 1, checkIfExists: true) { file -> (file.getSimpleName() =~ params.vcfchr_regex)[0][1] }
   .ifEmpty { error "Unable to find VCFs matching glob: ${params.vcf_glob}" }
   .tap { input_vcfs }
   .subscribe { println "Added chr${it[0]} (${it[1]}) to input_vcfs channels" }

Channel
   .fromFilePairs(params.vcf_glob+'.tbi', size: 1, checkIfExists: true) { file -> (file.getSimpleName() =~ params.vcfchr_regex)[0][1] }
   .ifEmpty { error "Unable to find VCF indices matching glob: ${params.vcf_glob}.tbi" }
   .tap { input_tbis }
   .subscribe { println "Added chr${it[0]} (${it[1]}) to input_tbis channels" }

bychrom_vcfs = input_vcfs.join(input_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
   .tap { bychrom_vcfs_smplstats }
   .tap { bychrom_vcfs_triostats }
   .tap { bychrom_vcfs_gtcheck }

//Set up the file channels for the ref and its various index components:
//Inspired by the IARC alignment-nf pipeline
//fai is generated by samtools faidx
ref = file(params.ref, checkIfExists: true)
ref_fai = file(params.ref+'.fai', checkIfExists: true)

//Set up the file channels for the trio PED, array VCF, and array sites TSV:
trioped = file(params.trio_ped, checkIfExists: true)
arrayvcf = file(params.array_vcf, checkIfExists: true)
arraytbi = file(params.array_vcf+'.tbi', checkIfExists: true)
arraysites = file(params.array_sites, checkIfExists: true)

//Defaults for cpus, memory, and time for each process:
//bcftools stats
params.stats_cpus = 1
params.stats_mem = 2
params.stats_timeout = '24h'
//bcftools +smpl-stats
params.smplstats_cpus = 1
params.smplstats_mem = 2
params.smplstats_timeout = '24h'
//bcftools +trio-stats
params.triostats_cpus = 1
params.triostats_mem = 2
params.triostats_timeout = '24h'
//bcftools gtcheck
params.gtcheck_cpus = 1
params.gtcheck_mem = 2
params.gtcheck_timeout = '24h'

process stats {
   tag "chr${chrom}"

   cpus params.stats_cpus
   memory { params.stats_mem.plus(task.attempt.minus(1).multiply(4))+' GB' }
   time { task.attempt == 2 ? '72h' : params.stats_timeout }
   errorStrategy { task.exitStatus in ([1]+(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(vcf), path(tbi) from bychrom_vcfs
   path ref
   path ref_fai

   output:
   tuple path("bcftools_stats_chr${chrom}_${params.run_suffix}.stderr"), path("bcftools_stats_chr${chrom}_${params.run_suffix}.stdout") into stats_logs

   shell:
   '''
   module load !{params.mod_bcftools}
   bcftools stats -I -F !{ref} !{params.filterstr} !{vcf} 2> bcftools_stats_chr!{chrom}_!{params.run_suffix}.stderr > bcftools_stats_chr!{chrom}_!{params.run_suffix}.stdout
   '''
}

process smplstats {
   tag "chr${chrom}"

   cpus params.smplstats_cpus
   memory { params.smplstats_mem.plus(task.attempt.minus(1).multiply(4))+' GB' }
   time { task.attempt == 2 ? '72h' : params.smplstats_timeout }
   errorStrategy { task.exitStatus in ([1]+(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(vcf), path(tbi) from bychrom_vcfs_smplstats
   path ref
   path ref_fai

   output:
   tuple path("bcftools_smplstats_chr${chrom}_${params.run_suffix}.stderr"), path("bcftools_smplstats_chr${chrom}_${params.run_suffix}.stdout") into smplstats_logs

   shell:
   '''
   module load !{params.mod_bcftools}
   bcftools view -Ov !{params.filterstr} !{vcf} | \
      bcftools +smpl-stats !{params.includestr} - 2> bcftools_smplstats_chr!{chrom}_!{params.run_suffix}.stderr > bcftools_smplstats_chr!{chrom}_!{params.run_suffix}.stdout
   '''
}

process triostats {
   tag "chr${chrom}"

   cpus params.triostats_cpus
   memory { params.triostats_mem.plus(task.attempt.minus(1).multiply(4))+' GB' }
   time { task.attempt == 2 ? '72h' : params.triostats_timeout }
   errorStrategy { task.exitStatus in ([1]+(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(vcf), path(tbi) from bychrom_vcfs_triostats
   path ref
   path ref_fai
   path trioped

   output:
   tuple path("bcftools_triostats_chr${chrom}_${params.run_suffix}.stderr"), path("bcftools_triostats_chr${chrom}_${params.run_suffix}.stdout") into triostats_logs

   shell:
   '''
   module load !{params.mod_bcftools}
   bcftools view -Ov !{params.filterstr} !{vcf} | \
      bcftools +trio-stats -p !{trioped} !{params.includestr} - 2> bcftools_triostats_chr!{chrom}_!{params.run_suffix}.stderr > bcftools_triostats_chr!{chrom}_!{params.run_suffix}.stdout
   '''
}

process gtcheck {
   tag "chr${chrom}"

   cpus params.gtcheck_cpus
   memory { params.gtcheck_mem.plus(task.attempt.minus(1).multiply(4))+' GB' }
   time { task.attempt == 2 ? '72h' : params.gtcheck_timeout }
   errorStrategy { task.exitStatus in ([1]+(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(vcf), path(tbi) from bychrom_vcfs_gtcheck
   path ref
   path ref_fai
   path arrayvcf
   path arraytbi
   path arraysites

   output:
   tuple path("bcftools_gtcheck_chr${chrom}_${params.run_suffix}.stderr"), path("bcftools_gtcheck_chr${chrom}_${params.run_suffix}.stdout") into gtcheck_logs

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   bcftools view -Oz -o query_subset_!{params.run_suffix}.vcf.gz -R !{arraysites} -T !{arraysites} !{params.filterstr} !{vcf}
   tabix -f query_subset_!{params.run_suffix}.vcf.gz
   bcftools gtcheck -u GT,GT -e 0 -g !{arrayvcf} query_subset_!{params.run_suffix}.vcf.gz 2> bcftools_gtcheck_chr!{chrom}_!{params.run_suffix}.stderr > bcftools_gtcheck_chr!{chrom}_!{params.run_suffix}.stdout
   '''
}

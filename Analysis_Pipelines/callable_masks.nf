#!/usr/bin/env nextflow
/* Pipeline to identify callable sites for each sample from BAMs            *
 * Core steps:                                                              *
 *  mosdepth equivalent to GATK CallableLoci on each BAM                    */

//Default paths, globs, and regexes:
//Glob for the per-individual BAMs and their indices:
params.bam_glob = "${projectDir}/BAMs/*_MD_IR_recal_filtered.ba{m,m.bai}"
//Regex for parsing the sample ID out of the BAM filename:
params.bam_regex = ~/^(.+)_MD_IR_recal_filtered$/
//Glob for the per-individual autosomal depth distributions from mosdepth
// generated by filter_VCF.nf:
params.depth_dist_glob = "${projectDir}/depthdists/*.mosdepth.global.dist.txt"
//Regex for parsing the sample ID out of the depth dist filename:
params.depth_regex = ~/^(.+)$/

//Set up the channel for all of the per-sample BAMs and their indices:
Channel
   .fromFilePairs(params.bam_glob, checkIfExists: true) { file -> (file.getSimpleName() =~ params.bam_regex)[0][1] }
   .ifEmpty { error "Unable to find BAMs matching glob: ${params.bam_glob}" }
   .tap { bams }
   .subscribe { println "Added ${it[0]} to bams channel" }
//And a channel for the mosdepth *.global.dist.txt files from filter_VCF.nf:
Channel
   .fromPath(params.depth_dist_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find depth distribution files matching glob: ${params.depth_dist_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.depth_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract sample ID from depth distribution filename: ${params.depth_regex}" }
   .tap { depth_dists }
   .subscribe { println "Added ${it[0]} (${it[1]}) to depth_dists channel" }

//Default parameter values:
//MAPQ threshold for calculating depth with mosdepth for GATK CallableLoci-like thresholding on DP:
params.mapq_thresh = "30"
//Lower and upper bounds for depth when making the per-sample masks:
//Set the *_threshold variables to blank ("") to use the quantiles instead
//It's allowed to use one threshold and one quantile
params.mindp_threshold = "10"
//params.mindp_quantile = "0.005"
params.mindp_quantile = ""
//params.maxdp_threshold = "80"
params.maxdp_threshold = ""
params.maxdp_quantile = "0.995"

//Defaults for cpus, memory, and time for each process:
//Constructing per-population mask BEDs
params.mosdepth_cpus = 1
params.mosdepth_mem = 2
params.mosdepth_timeout = '12h'

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

   output:
   tuple path("mosdepth_quantize_${sampleid}.stderr"), path("mosdepth_quantize_${sampleid}.stdout") into mosdepth_logs
   tuple val(sampleid), path("${sampleid}.mask.bed.gz") into mosdepth_masks
   tuple val(sampleid), path("${sampleid}.callable.bed.gz") into mosdepth_callable

   shell:
   '''
   module load !{params.mod_mosdepth}
   if [[ -z "!{params.mindp_quantile}" ]] && [[ -n "!{params.mindp_threshold}" ]]; then
      mindp=!{params.mindp_threshold}
   elif [[ -n "!{params.mindp_quantile}" ]]; then
      mindp=$(awk 'BEGIN{FS="\t";OFS=FS;mindpq=1-!{params.mindp_quantile};}$1=="total"&&$3<=mindpq{dpthresh=$2;}END{print dpthresh;}' !{depthdist})
   else
      echo "Error: Invalid specification of mindp_quantile or mindp_threshold, at least one must be set." > /dev/stderr
      exit 2
   fi
   if [[ -z "!{params.maxdp_quantile}" ]] && [[ -n "!{params.maxdp_threshold}" ]]; then
      maxdp=!{params.maxdp_threshold}
   elif [[ -n "!{params.maxdp_quantile}" ]]; then
      maxdp=$(awk 'BEGIN{FS="\t";OFS=FS;maxdpq=1-!{params.maxdp_quantile};}$1=="total"&&$3<=maxdpq{dpthresh=$2;}END{print dpthresh;}' !{depthdist})
   else
      echo "Error: Invalid specification of maxdp_quantile or maxdp_threshold, at least one must be set." > /dev/stderr
      exit 3
   fi
   MOSDEPTH_Q0=NO_DEPTH MOSDEPTH_Q1=LOW_DEPTH MOSDEPTH_Q2=CALLABLE MOSDEPTH_Q3=HIGH_DEPTH mosdepth -n --fast-mode --mapq !{params.mapq_thresh} --quantize 0:1:${mindp}:${maxdp}: !{sampleid}_quantized !{bambai[0]} 2> mosdepth_quantize_!{sampleid}.stderr > mosdepth_quantize_!{sampleid}.stdout
   gzip -dc !{sampleid}_quantized.quantized.bed.gz | \
      awk 'BEGIN{FS="\t";OFS=FS;}$4!="CALLABLE"{print $1, $2, $3;}' | \
      gzip -9 > !{sampleid}.mask.bed.gz
   gzip -dc !{sampleid}_quantized.quantized.bed.gz | \
      awk 'BEGIN{FS="\t";OFS=FS;}$4=="CALLABLE"{print $1, $2, $3;}' | \
      gzip -9 > !{sampleid}.callable.bed.gz
   '''
}

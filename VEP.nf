#!/usr/bin/env nextflow
/* Pipeline to run VEP on a whole-genome variant call set.    *
 * Core steps:                                                *
 *  Determine equal-density bands for each chromosome ->      *
 *  Run VEP on each equal-density band ->                     *
 *  Gather VEP results into per-chromosome outputs            */

nextflow.enable.dsl=2

//Default paths, globs, and regexes:
//Input VCFs:
params.vcf_glob = "${projectDir}/VCFs/*.vcf.gz"
params.vcf_regex = ~/_chr(\p{Alnum}+)_phased$/
//Band variant density for splitting:
params.band_num_variants = "250000"
//VEP parameters:
params.vep_cache = "/gpfs/gibbs/pi/tucci/pfr8/vep_r107_cache/"
params.vep_species = "homo_sapiens"
params.vep_assembly = "GRCh37"

//Defaults for cpus, base memory, memory ramp, memory unit, base timeout, timeout ramp, and timeout unit
//Calculate VCF splits
params.calcsplit_cpus = 1
params.calcsplit_mem = 1
params.calcsplit_mem_ramp = 4
params.calcsplit_mem_unit = 'GB'
params.calcsplit_timeout = 1
params.calcsplit_timeout_ramp = 4
params.calcsplit_timeout_unit = 'h'
//VEP
params.vep_cpus = 1
params.vep_mem = 2
params.vep_mem_ramp = 8
params.vep_mem_unit = 'GB'
params.vep_timeout = 1
params.vep_timeout_ramp = 4
params.vep_timeout_unit = 'h'
//Gather VEP VCFs
params.concatvcf_cpus = 1
params.concatvcf_mem = 2
params.concatvcf_mem_ramp = 8
params.concatvcf_mem_unit = 'GB'
params.concatvcf_timeout = 1
params.concatvcf_timeout_ramp = 4
params.concatvcf_timeout_unit = 'h'

process calcsplit {
   tag "chr${chrom}"

   cpus params.calcsplit_cpus
   memory { params.calcsplit_mem.plus(task.attempt.minus(1).multiply(params.calcsplit_mem_ramp))+' '+params.calcsplit_mem_unit }
   time { params.calcsplit_timeout.plus(task.attempt.minus(1).multiply(params.calcsplit_timeout_ramp))+params.calcsplit_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/bands", mode: 'copy', pattern: '*.bed'

   input:
   tuple val(chrom), path(vcf), path(tbi)

   output:
   tuple val(chrom), path("${params.run_name}_chr${chrom}_equaldensity_bands.bed")

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   bcftools query -f '%CHROM\\t%POS\\n' !{vcf} | \\
      !{projectDir}/HumanPopGenScripts/VEP/equalDensityBands.awk -v "bandsize=!{params.band_num_variants}" > !{params.run_name}_chr!{chrom}_equaldensity_bands.bed
   '''
}

process vep {
   tag "chr${chrom}b${band}"

   cpus params.vep_cpus
   memory { params.vep_mem.plus(task.attempt.minus(1).multiply(params.vep_mem_ramp))+' '+params.vep_mem_unit }
   time { params.vep_timeout.plus(task.attempt.minus(1).multiply(params.vep_timeout_ramp))+params.vep_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(chrom), val(band), path(bandbed), path(vcf), path(tbi)

   output:
   tuple val(chrom), path("chr${chrom}_${band}_VEP.vcf.gz"), path("chr${chrom}_${band}_VEP.vcf.gz.tbi")

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   module load !{params.mod_vep}
   #Compress and index the band BED:
   bgzip -c !{bandbed} > !{bandbed}.gz
   tabix -p bed !{bandbed}.gz
   #Extract the band from the full chromosome VCF:
   bcftools view -G -R !{bandbed}.gz -T !{bandbed}.gz -Oz -o chr!{chrom}_!{band}.vcf.gz !{vcf}
   tabix -f chr!{chrom}_!{band}.vcf.gz
   #Run VEP:
   cache_options="--cache --dir_cache !{params.vep_cache} --species !{params.vep_species} --assembly !{params.vep_assembly}"
   vepout_options="--vcf --regulatory --symbol --canonical --biotype --domains --no_stats --warning_file chr!{chrom}_!{band}_VEP_warnings.txt"
   vep ${cache_options} -i chr!{chrom}_!{band}.vcf.gz -o STDOUT ${vepout_options} | \\
      bgzip -c > chr!{chrom}_!{band}_VEP.vcf.gz
   tabix -f chr!{chrom}_!{band}_VEP.vcf.gz
   '''
}

process concatvcf {
   tag "chr${chrom}"

   cpus params.concatvcf_cpus
   memory { params.concatvcf_mem.plus(task.attempt.minus(1).multiply(params.concatvcf_mem_ramp))+' '+params.concatvcf_mem_unit }
   time { params.concatvcf_timeout.plus(task.attempt.minus(1).multiply(params.concatvcf_timeout_ramp))+params.concatvcf_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/VEP_VCFs", mode: 'copy', pattern: '*.vcf.g{z,z.tbi}'

   input:
   tuple val(chrom), path(vepvcfs)

   output:
   tuple val(chrom), path("${params.run_name}_chr${chrom}_VEP.vcf.gz"), path("${params.run_name}_chr${chrom}_VEP.vcf.gz.tbi")

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   #Stage the sharded VCFs:
   mkdir sharded_VCFs
   while IFS=$'\\t' read -a a;
      do
      ln -s ${a[0]} sharded_VCFs/${a[1]};
      ln -s ${a[2]} sharded_VCFs/${a[3]};
      echo "sharded_VCFs/${a[1]}";
   done < !{vepvcfs} | \\
      sort -k1,1V > chr!{chrom}_sharded_VEP_VCFs.fofn
   #Concatenate the sharded VCFs in order and index:
   bcftools concat -n -f chr!{chrom}_sharded_VEP_VCFs.fofn -Oz -o !{params.run_name}_chr!{chrom}_VEP.vcf.gz
   tabix -f !{params.run_name}_chr!{chrom}_VEP.vcf.gz
   '''
}

workflow {
   //Load the per-chromosome VCFs:
   Channel
      .fromPath(params.vcf_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find VCFs matching glob: ${params.vcf_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from VCF: ${params.vcf_regex}" }
      .tap { vcfs_only }
      .subscribe { println "Added ${it[1]} to vcfs_only channel" }
   Channel
      .fromPath(params.vcf_glob+'.tbi', checkIfExists: true)
      .ifEmpty { error "Unable to find VCF indices matching glob: ${params.vcf_glob}.tbi" }
      .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from VCF index: ${params.vcf_regex}" }
      .tap { tbis_only }
      .subscribe { println "Added ${it[1]} to tbis_only channel" }
   vcfs = vcfs_only
      .join(tbis_only, by: 0, failOnDuplicate: true, failOnMismatch: true)

   //Determine the equal-density bands for each chromosome:
   calcsplit(vcfs)

   //Perform VEP annotation of the variants:
   vep_inputs = calcsplit.out
      .splitText(by: 1, file: true, elem: 1)
      .map({ [ it[0], (it[1].getBaseName() =~ /[.]([0-9]+)$/)[0][1], it[1] ] })
      .combine(vcfs, by: 0)
   vep(vep_inputs)

   //Gather the VEP shards by chromosome
   concat_inputs = vep.out
      .collectFile() { [ "chr${it[0]}.txt", it[1].toAbsolutePath().toString()+'\t'+it[1].getName()+'\t'+it[2].toAbsolutePath().toString()+'\t'+it[2].getName()+'\n' ] }
      .map({ [ (it.getSimpleName() =~ /^chr(\p{Alnum}+)$/)[0][1], it ] })
   concatvcf(concat_inputs)
}

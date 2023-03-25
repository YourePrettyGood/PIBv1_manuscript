#!/usr/bin/env nextflow
/* Pipeline to perform projection of Sprime tracts onto phased genotypes    *
 * Core steps:                                                              *
 *  Project S' tracts onto haplotypes                                       */

//Default paths, globs, and regexes:
//Phased VCFs:
params.vcf_glob = "${projectDir}/phased_VCFs/*.vcf.gz"
//Regex for extracting chromosome from VCF filename:
params.vcf_regex = ~/_chr(\p{Alnum}+)$/

//Include/filter expression to apply to input VCFs:
params.input_filter_str = ""
//Sample ID file for exclusion:
params.samples_to_exclude = "<(echo)"

//Reference-related parameters for the pipeline:
params.ref_prefix = "/gpfs/gibbs/pi/tucci/pfr8/refs"
params.ref = "${params.ref_prefix}/1kGP/hs37d5/hs37d5.fa"
params.autosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
params.nuclearchroms = "${params.autosomes},X,Y"

//Set up the channels of per-chromosome phased VCFs and their indices:
Channel
   .fromPath(params.vcf_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find VCFs matching glob: ${params.vcf_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from VCF: ${params.vcf_regex}" }
   .tap { phased_vcfs }
   .subscribe { println "Added ${it[1]} to phased_vcfs channel" }

Channel
   .fromPath(params.vcf_glob+'.tbi', checkIfExists: true)
   .ifEmpty { error "Unable to find VCF indices matching glob: ${params.vcf_glob}.tbi" }
   .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from VCF index: ${params.vcf_regex}" }
   .tap { phased_tbis }
   .subscribe { println "Added ${it[1]} to phased_tbis channel" }

//Set up the channel of Sprime .score files (per-population and per-chromosome):
Channel
   .fromPath(params.sprime_score_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find Sprime score files matching glob: ${params.sprime_score_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.score_regex)[0][1], (a.getSimpleName() =~ params.score_regex)[0][2], a] }
   .ifEmpty { error "Regex failed to extract population name and chromosome ID from Sprime score file: ${params.score_regex}" }
   .tap { sprime_scores }
   .subscribe { println "Added ${it[2]} to sprime_scores channel" }

//Set up the channel of target populations for Sprime:
Channel
   .fromPath(params.target_pops_file, checkIfExists: true)
   .ifEmpty { error "Unable to find target populations file for Sprime: ${params.target_pop_file}" }
   .splitCsv(sep: "\t")
   .map { it[0] }
   .tap { sprime_pops_for_project }
   .tap { sprime_pops_for_cat }
   .subscribe { println "Added ${it} to sprime_pops_for_* channels" }

num_autosomes = params.autosomes.tokenize(',').size()

//Set up the file channels for the metadata files and the annotation GFF:
metadata = file(params.metadata_file, checkIfExists: true)
//And a file channel for the sample IDs to exclude:
excluded_samples = file(params.samples_to_exclude, checkIfExists: true)

//Default parameter values:
//Sample ID column name in metadata file:
params.id_colname = "Sample"
//Target group column name:
params.sprime_target_colname = "AnalysisGroup"
//Number of consecutive modern alleles allowed before splitting a tract:
//Slightly deceptively named
params.tract_max_gap = 5

//Defaults for cpus, memory, and time for each process:
//Sprime projection
params.sprimeproject_cpus = 1
params.sprimeproject_mem = 1
params.sprimeproject_timeout = '24h'
//Combine outputs across targets and chromosomes
params.catouts_cpus = 1
params.catouts_mem = 1
params.catouts_timeout = '2h'

//Preprocess the phased VCF channel to include the indices:
phased_vcfs
   .join(phased_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
   .filter({ params.autosomes.tokenize(',').contains(it[0]) })
   .set { phased_vcfs_for_project }

process sprime_project {
   tag "${pop} chr${chrom}"

   cpus params.sprimeproject_cpus
   memory { params.sprimeproject_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.sprimeproject_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

//   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.bed'

   input:
   tuple val(pop), val(chrom), path(score), path(vcf), path(tbi) from sprime_scores
      .join(sprime_pops_for_project.combine(phased_vcfs_for_project), by: [0,1], failOnDuplicate: true, failOnMismatch: true)
   path metadata
   path excluded_samples

   output:
   path("${pop}_chr${chrom}_Sprime_phased_tracts_perSample_maxgap${params.tract_max_gap}.bed") into sprime_project_BEDs

   shell:
   '''
   module load !{params.mod_bcftools}
   #Identify the samples in the target population:
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{pop}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_samples.tsv
   #Extract the set of Sprime sites for this chromosome:
   awk -v "c=!{chrom}" 'BEGIN{FS="\t";OFS=FS;}NR>1&&$1==c{print $1, $2;}' !{score} > !{pop}_chr!{chrom}_Sprime_sites.tsv
   #Project Sprime alleles onto each haplotype to identify tracts:
   bcftools query -f '%CHROM:%POS[\t%GT]\n' -T !{pop}_chr!{chrom}_Sprime_sites.tsv -H -S !{pop}_samples.tsv !{vcf} | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimePerSampleTracts.awk -v "spop=!{pop}" -v "allout=1" -v "phased=1" !{score} - | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBED.awk -v "phased=1" -v "max_gap=!{params.tract_max_gap}" | \
      sort -k1,1V -k2,2n -k3,3n > !{pop}_chr!{chrom}_Sprime_phased_tracts_perSample_maxgap!{params.tract_max_gap}.bed
   '''
}

process cat_outs {
   cpus params.catouts_cpus
   memory { params.catouts_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.catouts_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.t{sv,csv}.gz'

   input:
   path projections from sprime_project_BEDs.collectFile() { [ "projection_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path targetpops from sprime_pops_for_cat.collectFile(name: 'Sprime_target_populations.txt', newLine: true, sort: true)
   path metadata
   path excluded_samples

   output:
   path("${params.run_name}_Sprime_perChrom_perHaplotype_perPop_tract_lengths_maxgap${params.tract_max_gap}.tsv.gz") into catouts_outputs

   shell:
   '''
   #Symlink the many input files to avoid SLURM SBATCH script size limits:
   cat projection_paths.tsv | \
      while IFS=$'\t' read -a a;
         do
         ln -s ${a[2]} ${a[1]};
      done
   #Calculate total per-haplotype tract length and
   # concatenate files across populations:
   autosomes=!{params.autosomes}
   header=1
   while IFS=$'\t' read p;
      do
      for c in ${autosomes//,/$IFS};
         do
         !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=${p}" !{metadata} | \
            !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - | \
            !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBEDtoLengths.awk -v "header=${header}" -v "pop=${p}" -v "phased=1" -v "addpop=0" - ${p}_chr${c}_Sprime_phased_tracts_perSample_maxgap!{params.tract_max_gap}.bed
         header="";
      done;
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_perChrom_perHaplotype_perPop_tract_lengths_maxgap!{params.tract_max_gap}.tsv.gz
   unset header
   '''
}

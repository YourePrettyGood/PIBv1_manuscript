#!/usr/bin/env nextflow
/* Pipeline to identify and characterize introgression tracts from archaic  *
 *  hominins (e.g. Neandertals and Denisovans)                              *
 * Core steps:                                                              *
 *  (Reformat genetic/recombination rate map files for AS2 &&               *
 *  Extraction of target populations + outgroup from main VCF for AS2) ->   *
 *  ArchaicSeeker2 on each target population |->                            *
 *  MultiWaver2.1 on each target population for each archaic origin         *
 *  |-> match rate calculation                                              *
 *  |-> overlapping genes                                                   */

//Default paths, globs, and regexes:
//Jointly genotyped VCFs:
params.vcf_glob = "${projectDir}/perchrom_VCFs/*.vcf.gz"
//Regex for extracting chromosome from VCF filename:
params.vcf_regex = ~/_chr(\p{Alnum}+)$/
//Recombination rate maps:
params.recmap_glob = "${projectDir}/genetic_maps/genetic_map_GRCh37_chr*.txt"
//Regex for extracting chromosome from recombination rate map filename:
params.recmap_regex = ~/_chr([0-9]+)$/
//Archaic VCFs:
params.arcvcf_glob = "${projectDir}/archaic_VCFs/*.vcf.gz"
//PanTro FASTA glob:
params.pantro_glob = "${projectDir}/chr*.hg19.chimp.fa.gz"
//Regex for extracting chromosome from PanTro FASTA filename:
params.pantro_regex = ~/chr([0-9]+)/
//Ancestral allele FASTA glob:
params.anc_glob = "${projectDir}/homo_sapiens_ancestor_*.fa.gz"
//Regex for extracting chromosome from ancestral allele FASTA filename:
params.anc_regex = ~/^homo_sapiens_ancestor_([0-9]+)$/

//Include/filter expression to apply to input VCFs:
params.input_filter_str = ""
//Sample ID file for exclusion:
params.samples_to_exclude = "<(echo)"

//Reference-related parameters for the pipeline:
params.ref_prefix = "/gpfs/gibbs/pi/tucci/pfr8/refs"
params.ref = "${params.ref_prefix}/1kGP/hs37d5/hs37d5.fa"
params.autosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
params.nuclearchroms = "${params.autosomes},X,Y"

//Set up a value channel just of the types of archaic origins we want to run through MultiWaver2.1:
Channel
   .of("Neanderthal", "Denisova")
   .tap { archaics }
   .subscribe { println "Adding ${it} to archaics channel" }

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

//Set up the channel of recombination rate map files for ArchaicSeeker2:
Channel
   .fromPath(params.recmap_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find recombination rate map files matching glob: ${params.recmap_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.recmap_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from VCF: ${params.recmap_glob} =~ ${params.recmap_regex}" }
   .tap { recmaps }
   .subscribe { println "Added chr${it[0]} recombination rate map (${it[1]}) to recmaps channel" }

//Set up the channel of outgroup/Chimp allele FASTAs:
Channel
   .fromPath(params.pantro_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find outgroup/Chimp allele FASTAs matching glob: ${params.pantro_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.pantro_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from outgroup allele FASTA: ${params.pantro_glob} =~ ${params.pantro_regex}" }
   .tap { pantro }
   .subscribe { println "Added chr${it[0]} outgroup allele FASTA (${it[1]}) to pantro channel" }

//Set up the channel of ancestral allele FASTAs:
Channel
   .fromPath(params.anc_glob, checkIfExists: true)
   .ifEmpty { error "Unable to find ancestral allele FASTAs matching glob: ${params.anc_glob}" }
   .map { a -> [ (a.getSimpleName() =~ params.anc_regex)[0][1], a] }
   .ifEmpty { error "Regex failed to extract chromosome ID from ancestral allele FASTA: ${params.anc_glob} =~ ${params.anc_regex}" }
   .tap { ancestral }
   .subscribe { println "Added chr${it[0]} ancestral allele FASTA (${it[1]}) to ancestral channel" }

//Set up the channel of target populations for ArchaicSeeker2:
Channel
   .fromPath(params.target_pops_file, checkIfExists: true)
   .ifEmpty { error "Unable to find target populations file for ArchaicSeeker2: ${params.target_pop_file}" }
   .splitCsv(sep: "\t")
   .map { it[0] }
   .tap { as_pops }
   .tap { as_pops_for_cat }
   .subscribe { println "Added ${it} to as_pops channel" }

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
params.AS_outgroup_colname = "Region"
params.AS_outgroup = "Africa"
//Target group column name:
params.AS_target_colname = "AnalysisGroup"
//Archaic sample IDs:
//Options for NEA: AltaiNeandertal, Vindija33.19, Chagyrskaya-Phalanx
params.nea_id = "AltaiNeandertal"
params.den_id = "Denisova"
//AS2 model Newick tree:
//Note that this default tree is based on Altai Neanderthal, not Vindija or Chagyrskaya
params.AS2_model_tree = '((Outgroup:100,Modern:100):557.5,(Denisova:340,Neanderthal:300):237.5);'
//AS2 initial proportion of introgression:
//Default: 0.02
params.AS2_intro_proportion = "0.02"
//AS2 introgression time in generations:
//Default: 2000
params.AS2_intro_time = "2000"
//AS2 emission probability initial value:
//Default: 0.99
params.AS2_emission_prob = "0.99"
//Use all MW2 defaults, don't specify any optional arguments:
params.MW2_defaults = true
//MW2 lower bound for tract length:
//Default: 0
params.MW2_lower_bound = "0"
//MW2 significance threshold for the LRT:
//Default: 0.001
params.MW2_LRT_alpha = "0.001"
//MW2 number of bootstraps:
//Default: 100
params.MW2_n_bootstraps = "100"
//MW2 tolerance epsilon for convergence:
//Default: 1.0e-6
params.MW2_tolerance = "1.0e-6"
//MW2 minimum survival proportion:
//Default: 0.05
params.MW2_min_surv_prop = "0.05"
//MW2 maximum number of iterations:
//Default: 10000
params.MW2_max_iter = "10000"

//Defaults for cpus, memory, and time for each process:
//Reformatting of recombination rate map files:
//Memory for this one is specified in MiB, rest are in GiB
params.recmap_cpus = 1
params.recmap_mem = 512
params.recmap_timeout = '2h'
//VCF subsetting for ArchaicSeeker
params.ASsubset_cpus = 1
params.ASsubset_mem = 4
params.ASsubset_timeout = '24h'
//ArchaicSeeker2.0
params.AS2_cpus = 1
params.AS2_mem = 8
params.AS2_timeout = '24h'
//MultiWaver2.1
params.MW2_cpus = 1
params.MW2_mem = 8
params.MW2_timeout = '24h'
//Match rates
params.matchrate_cpus = 1
params.matchrate_mem = 2
params.matchrate_timeout = '2h'
//Gene lists
params.genes_cpus = 1
params.genes_mem = 5
params.genes_timeout = '24h'
//Combine outputs across targets and chromosomes
params.catouts_cpus = 1
params.catouts_mem = 1
params.catouts_timeout = '2h'

//Preprocess the per-chromosome VCF channel to include the indices:
perchrom_vcfs
   .join(perchrom_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
   .filter({ params.autosomes.tokenize(',').contains(it[0]) })
   .tap { perchrom_vcfs_OGsubset }
   .tap { perchrom_vcfs_for_matchrate }
   .set { perchrom_vcfs_subset }

//Do the same for the archaic VCF channel:
archaic_vcfs
   .join(archaic_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
   .filter({ params.autosomes.tokenize(',').contains(it[0]) })
   .tap { archaic_vcfs_for_matchrate }
   .set { archaic_vcfs_autosomes }

//Retain only the autosomes for the recombination rate maps:
recmaps
   .filter({ params.autosomes.tokenize(',').contains(it[0]) })
   .set { recmaps_autosomes }

//Same for the outgroup/Chimp FASTAs:
pantro
   .filter({ params.autosomes.tokenize(',').contains(it[0]) })
   .set { pantro_autosomes }

//Same for ancestral allele FASTAs:
ancestral
   .filter({ params.autosomes.tokenize(',').contains(it[0]) })
   .set { ancestral_autosomes }

process recmap {
   tag "chr${chrom}"

   cpus params.recmap_cpus
   memory { params.recmap_mem.plus(task.attempt.minus(1).multiply(512))+' MB' }
   time { task.attempt >= 2 ? '24h' : params.recmap_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(chrom), path(recmap) from recmaps_autosomes

   output:
   tuple val(chrom), path("recmap_${chrom}.txt") into recmaps_processed,recmaps_for_multiwaver

   shell:
   '''
   #Remove the header line and drop the first column:
   #This applies to the genetic map files provided for older versions of SHAPEIT
   tail -n+2 !{recmap} | cut -f2- > recmap_!{chrom}.txt
   '''
}

process vcf_subset {
   tag "${pop} chr${chrom}"

   cpus params.ASsubset_cpus
   memory { params.ASsubset_mem.plus(task.attempt.minus(1).multiply(8))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.ASsubset_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*_samples.tsv', saveAs: { 'chr'+chrom+'_'+it }

   input:
   tuple val(pop), val(chrom), path(input_vcf), path(input_tbi) from as_pops.combine(perchrom_vcfs_subset)
   path metadata
   path excluded_samples

   output:
   tuple path("bcftools_view_selectpops_AS2_${pop}_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_AS2_${pop}_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_AS2_${pop}_chr${chrom}.stdout") into AS_vcf_subset_logs
   tuple val(pop), path("${pop}_chr${chrom}_nomissinggenos.vcf.gz"), path("${pop}_chr${chrom}_nomissinggenos.vcf.gz.tbi") into AS_perchrom_vcfs
   tuple val(pop), path("${pop}_samples.tsv") into AS_pop_samples

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.AS_target_colname}" -v "select=!{pop}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_samples.tsv
   #We have to separate the sample selection from the genotype missingness filter
   # because of a bug in bcftools view when -S and -g are used in tandem:
   #Way too many sites get filtered if you use them together, possibly more than
   # if -g was used alone.
   bcftools view -S !{pop}_samples.tsv -Ou !{input_vcf} 2> bcftools_view_selectpops_AS2_!{pop}_chr!{chrom}.stderr | \
   bcftools view -i '!{params.input_filter_str}' -g ^miss -Oz -o !{pop}_chr!{chrom}_nomissinggenos.vcf.gz 2> bcftools_view_nomissinggenos_AS2_!{pop}_chr!{chrom}.stderr > bcftools_view_nomissinggenos_AS2_!{pop}_chr!{chrom}.stdout
   tabix -f !{pop}_chr!{chrom}_nomissinggenos.vcf.gz
   '''
}

process OG_vcf_subset {
   tag "chr${chrom}"

   cpus params.ASsubset_cpus
   memory { params.ASsubset_mem.plus(task.attempt.minus(1).multiply(8))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.ASsubset_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*_samples.tsv', saveAs: { 'chr'+chrom+'_'+it }

   input:
   tuple val(chrom), path(input_vcf), path(input_tbi) from perchrom_vcfs_OGsubset
   path metadata
   path excluded_samples

   output:
   tuple path("bcftools_view_selectpops_AS2_OG_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_AS2_OG_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_AS2_OG_chr${chrom}.stdout") into AS_vcf_OGsubset_logs
   tuple path("OG_chr${chrom}_nomissinggenos.vcf.gz"), path("OG_chr${chrom}_nomissinggenos.vcf.gz.tbi") into AS_OG_perchrom_vcfs
   path("OG_samples.tsv") into AS_OG_samples

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.AS_outgroup_colname}" -v "select=!{params.AS_outgroup}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > OG_samples.tsv
   #We have to separate the sample selection from the genotype missingness filter
   # because of a bug in bcftools view when -S and -g are used in tandem:
   #Way too many sites get filtered if you use them together, possibly more than
   # if -g was used alone.
   bcftools view -S OG_samples.tsv -Ou !{input_vcf} 2> bcftools_view_selectpops_AS2_OG_chr!{chrom}.stderr | \
   bcftools view -i '!{params.input_filter_str}' -g ^miss -Oz -o OG_chr!{chrom}_nomissinggenos.vcf.gz 2> bcftools_view_nomissinggenos_AS2_OG_chr!{chrom}.stderr > bcftools_view_nomissinggenos_AS2_OG_chr!{chrom}.stdout
   tabix -f OG_chr!{chrom}_nomissinggenos.vcf.gz
   '''
}

process arc_vcf_subset {
   tag "chr${chrom}"

   cpus params.ASsubset_cpus
   memory { params.ASsubset_mem.plus(task.attempt.minus(1).multiply(8))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.ASsubset_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(input_vcf), path(input_tbi) from archaic_vcfs_autosomes

   output:
   tuple path("bcftools_split_AS2_arc_chr${chrom}.stderr"), path("bcftools_split_AS2_arc_chr${chrom}.stdout") into AS_vcf_arcsubset_logs
   tuple path("NEA_chr${chrom}.vcf.gz"), path("NEA_chr${chrom}.vcf.gz.tbi") into AS_NEA_perchrom_vcfs
   tuple path("DEN_chr${chrom}.vcf.gz"), path("DEN_chr${chrom}.vcf.gz.tbi") into AS_DEN_perchrom_vcfs

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   printf "!{params.nea_id}\t-\tNEA_chr!{chrom}.vcf.gz\n" > archaic_samples.tsv
   printf "!{params.den_id}\t-\tDEN_chr!{chrom}.vcf.gz\n" >> archaic_samples.tsv
   bcftools +split -S archaic_samples.tsv -Oz -o . !{input_vcf} 2> bcftools_split_AS2_arc_chr!{chrom}.stderr > bcftools_split_AS2_arc_chr!{chrom}.stdout
   tabix -f NEA_chr!{chrom}.vcf.gz
   tabix -f DEN_chr!{chrom}.vcf.gz
   '''
}

process AS {
   tag "${pop}"

   cpus params.AS2_cpus
   memory { params.AS2_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.AS2_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/ArchaicSeeker", mode: 'copy', pattern: '*.{seg,sum}'

   input:
   tuple val(pop), path(targetvcfs), path(othervcfs), path(recmaps), path(pantro), path(anc), path(OG), path(popsamples) from AS_perchrom_vcfs
      .collectFile() { [ "${it[0]}.txt", it[1].toAbsolutePath().toString()+'\t'+it[1].getName()+'\t'+it[2].toAbsolutePath().toString()+'\t'+it[2].getName()+'\n' ] }
      .map({ [ it.getSimpleName(), it ] })
      .combine(AS_OG_perchrom_vcfs
         .mix(AS_NEA_perchrom_vcfs)
         .mix(AS_DEN_perchrom_vcfs)
         .collectFile() { [ "other_vcfs.txt", it[0].toAbsolutePath().toString()+'\t'+it[0].getName()+'\t'+it[1].toAbsolutePath().toString()+'\t'+it[1].getName()+'\n' ] })
      .combine(recmaps_processed
         .collectFile() { [ "recmaps.txt", it[1].toAbsolutePath().toString()+'\t'+it[1].getName()+'\t'+it[0]+'\n' ] })
      .combine(pantro_autosomes
         .collectFile() { [ "pantro.txt", it[1].toAbsolutePath().toString()+'\t'+it[1].getName()+'\t'+it[0]+'\n' ] })
      .combine(ancestral_autosomes
         .collectFile() { [ "anc.txt", it[1].toAbsolutePath().toString()+'\t'+it[1].getName()+'\t'+it[0]+'\n' ] })
      .combine(AS_OG_samples.first())
      .join(AS_pop_samples
         .groupTuple(by: 0, size: num_autosomes)
         .map({ [it[0], it[1].take(1)] }), by: 0, failOnDuplicate: true, failOnMismatch: true)
   path metadata

   output:
   tuple path("AS2_${pop}.stderr"), path("AS2_${pop}.stdout") into AS_logs
   tuple val(pop), path("AS2_${pop}.seg"), path("AS2_${pop}.sum") into AS2_outputs,AS2_outputs_for_matchrate,AS2_outputs_for_genes

   shell:
   AS2_params = "-a ${params.AS2_intro_proportion} -T ${params.AS2_intro_time} -e ${params.AS2_emission_prob}"
   '''
   module load !{params.mod_ArchaicSeeker}
   #Stage the target pop VCFs:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[0]} ${a[1]};
      ln -s ${a[2]} ${a[3]};
   done < !{targetvcfs}
   #Stage the other VCFs:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[0]} ${a[1]};
      ln -s ${a[2]} ${a[3]};
   done < !{othervcfs}
   #Stage the rec rate maps and FASTA files:
   cat !{recmaps} !{pantro} !{anc} | \
      while IFS=$'\t' read -a a;
         do
         ln -s ${a[0]} ${a[1]};
      done
   #Compose the VCF parameters file:
   printf "vcf\n" > vcf.par
   find . -name "!{pop}_chr*.vcf.gz" -print | \
      sort -k1,1V >> vcf.par
   find . -name "OG_chr*.vcf.gz" -print | \
      sort -k1,1V >> vcf.par
   find . -name "NEA_chr*.vcf.gz" -print | \
      sort -k1,1V >> vcf.par
   find . -name "DEN_chr*.vcf.gz" -print | \
      sort -k1,1V >> vcf.par
   #Compose the outgroup FASTA parameters file:
   printf "outgroup\tcontig\n" > outgroup.par
   sort -k3,3V !{pantro} | cut -f2,3 >> outgroup.par
   #Compose the ancestral allele FASTA parameters file:
   printf "ancestor\tcontig\n" > anc.par
   sort -k3,3V !{anc} | cut -f2,3 >> anc.par
   #Compose the recombination rate map parameters file:
   printf "remap\tcontig\n" > remap.par;
   sort -k3,3V !{recmaps} | cut -f2,3 >> remap.par
   #Compose the population labels parameters file:
   printf "ID\tPop\tArchaicSeekerPop\n" > pop.par
   printf "!{params.nea_id}\tNeanderthal\tArchaic\n" >> pop.par
   printf "!{params.den_id}\tDenisova\tArchaic\n" >> pop.par
   awk 'BEGIN{FS="\t";OFS=FS;}FNR==NR{print $1, "Outgroup", "African";}FNR<NR{print $1, "Modern", "Test";}' !{OG} !{popsamples} >> pop.par
   #Create the ArchaicSeeker2 model Newick tree file:
   printf "!{params.AS2_model_tree}\n" > model.txt
   #Now run ArchaicSeeker2:
   ArchaicSeeker2 !{AS2_params} -v vcf.par -r remap.par -p pop.par -X outgroup.par -A anc.par -m model.txt -o AS2_!{pop} 2> AS2_!{pop}.stderr > AS2_!{pop}.stdout
   '''
}

process MW {
   tag "${pop} ${arc}"

   cpus params.MW2_cpus
   memory { params.MW2_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.MW2_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/ArchaicSeeker", mode: 'copy', pattern: '*.{seg,sum}'

   input:
   tuple val(pop), path(seg), path(sum), path(recmaps) from AS2_outputs
      .combine(recmaps_for_multiwaver
         .collectFile() { [ "recmaps.txt", it[1].toAbsolutePath().toString()+'\tgenetic_map_chr'+it[0]+'_combined_b37.txt\n' ] })
   path metadata
   each arc from archaics

   output:
   tuple path("MultiWaver2.1_${pop}_${arc}.stderr"), path("MultiWaver2.1_${pop}_${arc}.stdout") into MW2_logs
   tuple val(pop), path("AS2_${pop}_${arc}.seg"), path("MW2_${pop}_${arc}.bootstrap"), path("MW2_${pop}_${arc}.sum") into MW2_outputs

   shell:
   MW2_params = params.MW2_defaults ? "" : "-l ${params.MW2_lower_bound} -a ${params.MW2_LRT_alpha} -b ${params.MW2_n_bootstraps} -e ${params.MW2_tolerance} -p ${params.MW2_min_surv_prop} -m ${params.MW2_max_iter}"
   '''
   module load !{params.mod_ArchaicSeeker}
   #Stage the genetic/recombination rate map files with appropriate names:
   mkdir genmap
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[0]} genmap/${a[1]};
   done < !{recmaps}
   #Convert ArchaicSeeker2 tracts into MultiWaver2.1 input format:
   getAS2Seg genmap !{seg} AS2_!{pop}_!{arc} !{arc} 2> getAS2Seg_!{pop}_!{arc}.stderr > getAS2Seg_!{pop}_!{arc}.stdout
   #Run MultiWaver2.1 on the ArchaicSeeker2 tracts:
   MultiWaver2.1 !{MW2_params} -i AS2_!{pop}_!{arc}.seg -t !{task.cpus} -o MW2_!{pop}_!{arc} 2> MultiWaver2.1_!{pop}_!{arc}.stderr > MultiWaver2.1_!{pop}_!{arc}.stdout
   '''
}

process matchrate {
   tag "${pop} chr${chrom}"

   cpus params.matchrate_cpus
   memory { params.matchrate_mem.plus(task.attempt.minus(1).multiply(4))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.matchrate_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(pop), path(seg), val(chrom), path(modernvcf), path(moderntbi), path(arcvcf), path(arctbi) from AS2_outputs_for_matchrate
                                                                                                               .map({ [it[0], it[1]] })
                                                                                                               .combine(perchrom_vcfs_for_matchrate
                                                                                                                           .join(archaic_vcfs_for_matchrate, by: 0, failOnDuplicate: true, failOnMismatch: true))
   path metadata
   path excluded_samples
   path archaic_metadata

   output:
   path("AS2_${pop}_chr${chrom}_match_rates.tsv.gz") into matchrates
   path("AS2_${pop}_chr${chrom}_matches.tsv.gz") into matches

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   #Get the list of sample IDs for this population:
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.AS_target_colname}" -v "select=!{pop}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_samples.tsv
   #Convert AS2 .seg format to something resembling Sprime .score format:
   bcftools query -S !{pop}_samples.tsv -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' !{modernvcf} | \
      !{projectDir}/HumanPopGenScripts/ArchaicSeeker/segToSprimeScore.awk -v "chrom=!{chrom}" <(tail -n+2 !{seg} | sort -k2,2V -k3,3n -k4,4n -k1,1V | cat <(head -n1 !{seg}) -) - > AS2_!{pop}_chr!{chrom}.score
   #Extract the list of sites from this .score file:
   fgrep -v "CHROM" AS2_!{pop}_chr!{chrom}.score | cut -f1,2 | uniq | sort -k1,1V -k2,2n > !{pop}_chr!{chrom}_sites.tsv
   #Now identify matches to the archaics in the usual Sprime way:
   bcftools query -T !{pop}_chr!{chrom}_sites.tsv -H -f '%CHROM:%POS[\t%TGT]\n' !{arcvcf} | \
      !{projectDir}/HumanPopGenScripts/Sprime/archaicMatchSprime.awk -v "spop=!{pop}" -v "gtfmt=query" !{archaic_metadata} - AS2_!{pop}_chr!{chrom}.score | \
      gzip -9 > AS2_!{pop}_chr!{chrom}_matches.tsv.gz
   #And summarize each AS2 tract with match rates:
   gzip -dc AS2_!{pop}_chr!{chrom}_matches.tsv.gz | \
      sed 's/ASSTATE/SCORE/' | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeArchaicMatchRate.awk | \
      sed 's/SCORE/ASSTATE/' | \
      gzip -9 > AS2_!{pop}_chr!{chrom}_match_rates.tsv.gz
   '''
}

process genes {
   tag "${pop}"

   cpus params.genes_cpus
   memory { params.genes_mem.plus(task.attempt.minus(1).multiply(8))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.genes_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/ArchaicSeeker", mode: 'copy', pattern: '*.tcsv.gz'

   input:
   tuple val(pop), path(seg) from AS2_outputs_for_genes
                                     .map({ [it[0], it[1]] })
   path annotation_gff

   output:
   tuple path("${pop}_AS2_gene_lists.tcsv.gz"), path("${pop}_AS2_gene_name_lists.tcsv.gz") into gene_lists

   shell:
   '''
   module load !{params.mod_bedtools}
   #Identify genes overlapping AS2 tracts:
   bedtools intersect \
      -a <(tail -n+2 !{seg} | sort -k2,2V -k3,3n -k4,4n -k1,1V | cat <(head -n1 !{seg}) - | \
         !{projectDir}/HumanPopGenScripts/ArchaicSeeker/AS2segToBED.awk -v "pop=!{pop}" | \
         sort -k1,1 -k2,2n -k3,3n) \
      -b <(gzip -dc !{annotation_gff} | \
         awk 'BEGIN{FS="\t";OFS=FS;}/^#/{print;}!/^#/{sub("chr", "", $1); print $0;}') \
      -wao | \
      tee >(!{projectDir}/HumanPopGenScripts/Sprime/SprimeTractGeneList.awk -v "trim=1" -v "tag=gene_name" | \
         sort -k1,1V -k2,2n -k3,3n | \
         gzip -9 > !{pop}_AS2_gene_name_lists.tcsv.gz) | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractGeneList.awk -v "trim=1" -v "tag=ID" | \
      sort -k1,1V -k2,2n -k3,3n | \
      gzip -9 > !{pop}_AS2_gene_lists.tcsv.gz
   '''
}

process cat_outs {
   cpus params.catouts_cpus
   memory { params.catouts_mem.plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.catouts_timeout }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/ArchaicSeeker", mode: 'copy', pattern: '*.t{sv,csv}.gz'

   input:
   path matchrates from matchrates.collectFile() { [ "matchrate_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path matches from matches.collectFile() { [ "match_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path genelists from gene_lists.flatten().collectFile() { [ "genelist_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   path targetpops from as_pops_for_cat.collectFile(name: 'AS2_target_populations.txt', newLine: true, sort: true)
   path metadata
   path excluded_samples

   output:
   tuple path("${params.run_name}_perPop_AS2_autosomal_match_rates.tsv.gz"), path("${params.run_name}_perPop_AS2_autosomal_matches.tsv.gz"), path("${params.run_name}_AS2_gene_name_lists.tcsv.gz"), path("${params.run_name}_AS2_gene_lists.tcsv.gz") into catouts_outputs

   shell:
   '''
   #Symlink the many input files to avoid SLURM SBATCH script size limits:
   cat matchrate_paths.tsv match_paths.tsv genelist_paths.tsv | \
      while IFS=$'\t' read -a a;
         do
         ln -s ${a[2]} ${a[1]};
      done
   #Concatenate the match rate files across populations and chromosomes:
   autosomes=!{params.autosomes}
   header=1
   while IFS=$'\t' read p;
      do
      for c in ${autosomes//,/$IFS};
         do
         gzip -dc AS2_${p}_chr${c}_match_rates.tsv.gz | \
            awk -v "header=${header}" 'BEGIN{FS="\t";OFS=FS;}NR==1&&header{print;}NR>1{print;}'
         header=0;
      done
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_perPop_AS2_autosomal_match_rates.tsv.gz
   unset header
   #Concatenate the match files across populations and chromosomes:
   header=1
   while IFS=$'\t' read p;
      do
      for c in ${autosomes//,/$IFS};
         do
         gzip -dc AS2_${p}_chr${c}_matches.tsv.gz | \
            awk -v "header=${header}" 'BEGIN{FS="\t";OFS=FS;}NR==1&&header{print;}NR>1{print;}'
         header=0;
      done;
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_perPop_AS2_autosomal_matches.tsv.gz
   unset header
   #Concatenate the gene name lists across populations:
   header=1
   while IFS=$'\t' read p;
      do
      if [[ "${header}" == "1" ]]; then
         printf "Chromosome\tStart\tEnd\tTractID\tOverlappingGenes\n"
         header=""
      fi
      gzip -dc ${p}_AS2_gene_name_lists.tcsv.gz
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_AS2_gene_name_lists.tcsv.gz
   unset header
   #Concatenate the gene ID lists across populations:
   header=1
   while IFS=$'\t' read p;
      do
      if [[ "${header}" == "1" ]]; then
         printf "Chromosome\tStart\tEnd\tTractID\tOverlappingGenes\n"
         header=""
      fi
      gzip -dc ${p}_AS2_gene_lists.tcsv.gz
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_AS2_gene_lists.tcsv.gz
   unset header
   '''
}

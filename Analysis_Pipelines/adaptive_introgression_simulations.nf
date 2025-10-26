#!/usr/bin/env nextflow
/* Pipeline to evaluate performance of adaptive introgression detection     *
 *  on neutral simulations using the PapuansOutOfAfrica_10J19 model         *
 *  hominins (e.g. Neandertals and Denisovans)                              *
 * Core steps:                                                              *
 *  Simulation for each autosome using stdpopsim/msprime ->                 *
 *  Conversion of TreeSequence to VCF ->                                    *
 *  Splitting into archaics-only VCF and moderns-only VCF -1>               *
 *  Extraction of target pop + outgroup for Sprime ->                       *
 *  Concatenation of per-autosome VCFs for Sprime ->                        *
 *  Sprime for target pop -2>                                               *
 *  2> (Calculate match rates -3> && Project S' tracts onto individuals)    *
 *  1,2,3> Extract archaic-matching S' variant sites and cluster by r^2 ->  *
 *  Estimate tract frequencies                                              *
 *  1> (Calculate windowed PBS && Calculate XP-EHH) ->                      *
 *  Combine outputs across autosomes and pops                               */

nextflow.enable.dsl=2

//Default paths, globs, and regexes:
//Recombination rate maps:
params.recmap_glob = "${projectDir}/"
//Regex for extracting chromosome from recombination rate map filename:
params.recmap_regex = ~/chr([0-9XY_par]+)/

//Reference-related parameters for the pipeline:
params.autosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
params.nuclearchroms = "${params.autosomes},X,Y"

//Simulation parameters for stdpopsim:
params.sim_seed = 42
params.sim_species = "HomSap"
params.sim_gmap = "HapMapII_GRCh37"
params.sim_demo_model = "PapuansOutOfAfrica_10J19"
params.sim_modern_groups = "AFR,CEU,CHB,OCN"
params.sim_arc_groups = "NeaA,DenA"
params.n_AFR = 101
params.n_EUR = 150
params.n_EAS = 237
params.n_OCN = 100
params.n_DEN = 1
params.n_NEA = 1
params.sim_sample_string = "YRI:${params.n_AFR} CEU:${params.n_EUR} CHB:${params.n_EAS} Papuan:${params.n_OCN} DenA:${params.n_DEN} NeaA:${params.n_NEA}"

//Default parameter values:
//Sample ID column name in metadata file:
params.id_colname = "SampleID"
//Outgroup name to use:
params.sprime_outgroup_colname = "Population"
params.sprime_outgroup = "YRI"
//Outgroup and target subsets to make:
params.pops_to_subset = "YRI,Papuan"
//Subset sizes:
params.pop_subset_sizes = "21,25"
//Region names to assign to the simulated pops (including subsets):
params.subset_region_names = "AFR,OCN"
//Target group column name:
params.sprime_target_colname = "Population"
//Column name to use for r^2 calculation:
params.r2_colname = "Region"

//Archaic hominin groups represented in archaic VCFs:
params.archaic_groups = "Neandertal,Denisovan"
//Archaic hominin samples to consider individually in archaic VCFs:
params.archaic_indivs = ""
//Short names for these archaic hominin samples:
params.archaic_shortnames = ""

//Tract origins to assess and criteria to use for core haplotypes:
params.origin_criteria = "PFR"
params.tract_origin_list = "Neandertal,Denisovan,Ambiguous"

//Internal gap length in phased projection allowed to be spanned:
//This gap length is a count of Sprime sites at which the haplotype
// in question does not match the Sprime-identified putative archaic
// allele.
//Strongly recommend setting to 0, as this matches haplotype
// predictions from other methods better.
params.tract_max_gap = 0

//Whether r^2 should be calculated based on genotypes or haplotypes:
params.phased = false
if (params.phased) {
   params.r2_type = "hap"
} else {
   params.r2_type = "geno"
}
//Minimum r^2 to report:
params.min_r2 = "0.3"

//PBS window size (in # variants):
params.window_size = 20
//PBS window step (in # variants):
params.window_step = 5
//Extra parameter to pbs_cli.py if we want to clip negative PBS scores to 0:
//params.clip_pbs = "--clip_neg_pbs"
params.clip_pbs = ""
//XP-EHH genetic map column positions (0-based):
//Physical position column:
// For plink.chr*.GRCh37.map, this is 3
params.gmap_pos_col = 3
//Genetic map position column:
// For plink.chr*.GRCh37.map, this is 2
params.gmap_cm_col = 2

//Defaults for cpus, memory, and time for each process:
//Run simulations and generate VCFs:
params.run_sim_cpus = 1
params.run_sim_mem = 8
params.run_sim_mem_increment = 16
params.run_sim_timeout = '4h'
params.run_sim_retry_timeout = '24h'
params.run_sim_queue = 'day'
params.run_sim_retry_queue = 'day'
//VCF subsetting for Sprime
params.vcfsubset_cpus = 1
params.vcfsubset_mem = 4
params.vcfsubset_mem_increment = 8
params.vcfsubset_timeout = '2h'
params.vcfsubset_retry_timeout = '24h'
params.vcfsubset_queue = 'day'
params.vcfsubset_retry_queue = 'day'
//VCF concatenating for Sprime
params.concatvcf_cpus = 1
params.concatvcf_mem = 4
params.concatvcf_mem_increment = 16
params.concatvcf_timeout = '2h'
params.concatvcf_retry_timeout = '48h'
params.concatvcf_queue = 'day'
params.concatvcf_retry_queue = 'week'
//Sprime
params.sprime_cpus = 1
params.sprime_mem = 8
params.sprime_mem_increment = 16
params.sprime_timeout = '2h'
params.sprime_retry_timeout = '24h'
params.sprime_queue = 'day'
params.sprime_retry_queue = 'day'
//Sprime match rates
params.matchrate_cpus = 1
params.matchrate_mem = 4
params.matchrate_mem_increment = 16
params.matchrate_timeout = '4h'
params.matchrate_retry_timeout = '24h'
params.matchrate_queue = 'day'
params.matchrate_retry_queue = 'day'
//Sprime projection
params.projection_cpus = 1
params.projection_mem = 1
params.projection_mem_increment = 16
params.projection_timeout = '4h'
params.projection_retry_timeout = '24h'
params.projection_queue = 'day'
params.projection_retry_queue = 'day'
//Identify Sprime "tag" sites
params.tagsites_cpus = 1
params.tagsites_mem = 1
params.tagsites_mem_increment = 4
params.tagsites_timeout = '3h'
params.tagsites_retry_timeout = '24h'
params.tagsites_queue = 'day'
params.tagsites_retry_queue = 'day'
//Identify core haplotypes based on LD connected components
params.calcld_cpus = 1
params.calcld_mem = 4
params.calcld_mem_increment = 4
params.calcld_timeout = '6h'
params.calcld_retry_timeout = '24h'
params.calcld_queue = 'day'
params.calcld_retry_queue = 'day'
//Concatenate core haplotypes across chromosomes and tract origins
params.catcorehaps_cpus = 1
params.catcorehaps_mem = 1
params.catcorehaps_mem_increment = 4
params.catcorehaps_timeout = '3h'
params.catcorehaps_retry_timeout = '24h'
params.catcorehaps_queue = 'day'
params.catcorehaps_retry_queue = 'day'
//Annotate .score/_matches.tsv file with core haplotypes
params.corescore_cpus = 1
params.corescore_mem = 4
params.corescore_mem_increment = 4
params.corescore_timeout = '3h'
params.corescore_retry_timeout = '24h'
params.corescore_queue = 'day'
params.corescore_retry_queue = 'day'
//Sprime core haplotype frequencies
params.coretf_cpus = 1
params.coretf_mem = 1
params.coretf_mem_increment = 4
params.coretf_timeout = '12h'
params.coretf_retry_timeout = '24h'
params.coretf_queue = 'day'
params.coretf_retry_queue = 'day'
//Construct map from core haplotype to Sprime tract
params.tractmap_cpus = 1
params.tractmap_mem = 1
params.tractmap_mem_increment = 4
params.tractmap_timeout = '1h'
params.tractmap_retry_timeout = '24h'
params.tractmap_queue = 'day'
params.tractmap_retry_queue = 'day'
//PBS
params.pbs_cpus = 1
params.pbs_mem = 24
params.pbs_mem_increment = 16
params.pbs_timeout = '2h'
params.pbs_retry_timeout = '24h'
params.pbs_queue = 'day'
params.pbs_retry_queue = 'day'
//XP-EHH
params.xpehh_cpus = 1
params.xpehh_mem = 32
params.xpehh_mem_increment = 16
params.xpehh_timeout = '2h'
params.xpehh_retry_timeout = '24h'
params.xpehh_queue = 'day'
params.xpehh_retry_queue = 'day'
//Combine outputs across targets and chromosomes
params.catouts_cpus = 1
params.catouts_mem = 1
params.catouts_mem_increment = 16
params.catouts_timeout = '2h'
params.catouts_retry_timeout = '24h'
params.catouts_queue = 'day'
params.catouts_retry_queue = 'day'

process run_sim {
   tag "${chrom}"

   cpus params.run_sim_cpus
   memory { params.run_sim_mem.plus(task.attempt.minus(1).multiply(params.run_sim_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.run_sim_retry_timeout : params.run_sim_timeout }
   queue { task.attempt >= 2 ? params.run_sim_retry_queue : params.run_sim_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/sims", mode: 'copy', pattern: '*.{ts,vcf.gz,vcf.gz.tbi,tsv}'

   input:
   val(chrom)
   val(num_autosomes)
   val(autosome_map)

   output:
   tuple path("stdpopsim_${params.run_name}_${chrom}.stderr"), path("stdpopsim_${params.run_name}_${chrom}.stdout"), path("tskit_vcf_${params.run_name}_${chrom}.stderr"), path("bgzip_${params.run_name}_${chrom}.stderr"), path(""), path(""), emit: logs
   tuple val(chrom), path("${params.run_name}_${chrom}.ts"), emit: ts
   tuple val(chrom), path("${params.run_name}_${chrom}_modern.vcf.gz"), path("${params.run_name}_${chrom}_modern.vcf.gz.tbi"), emit: modern_vcf
   tuple val(chrom), path("${params.run_name}_${chrom}_archaic.vcf.gz"), path("${params.run_name}_${chrom}_archaic.vcf.gz.tbi"), emit: archaic_vcf
   tuple val(chrom), path("${params.run_name}_${chrom}_sample_metadata.tsv"), emit: metadata

   shell:
   chrid = autosome_map[chrom]
   '''
   module load !{params.mod_miniconda}
   conda activate adintsims
   ((seed=!{params.sim_seed}*!{num_autosomes}+!{chrom}))
   #Run the simulation using msprime:
   stdpopsim !{params.sim_species} -s ${seed} -c !{chrid} -g !{params.sim_gmap} -d !{params.sim_demo_model} -o !{params.run_name}_!{chrom}.ts !{params.sim_sample_string} 2> stdpopsim_!{params.run_name}_!{chrom}.stderr > stdpopsim_!{params.run_name}_!{chrom}.stdout
   #Convert TreeSequence to VCF:
   tskit vcf -c !{chrid} !{params.run_name}_!{chrom}.ts 2> tskit_vcf_!{params.run_name}_!{chrom}.stderr | \\
      bgzip 2> bgzip_!{params.run_name}_!{chrom}.stderr > !{params.run_name}_!{chrom}.vcf.gz
   #Prepare the metadata file:
   !{projectDir}/HumanPopGenScripts/Simulations/ts_to_metadata.awk <(tskit populations !{params.run_name}_!{chrom}.ts) <(tskit nodes !{params.run_name}_!{chrom}.ts) | \\
      !{projectDir}/HumanPopGenScripts/Simulations/addRegionSplitMetadata.awk -v "subsets=!{params.pops_to_subset}" -v "subsetsizes=!{params.pop_subset_sizes}" -v "regions=!{params.subset_region_names}" > !{params.run_name}_!{chrom}_sample_metadata.tsv
   #Partition the VCF into modern samples only and archaic samples only:
   !{projectDir}/HumanPopGenScripts/Simulations/metadata_to_split_map.awk -v "prefix=!{params.run_name}_!{chrom}" -v "modern=!{params.sim_modern_groups}" -v "arc=!{params.sim_arc_groups}" !{params.run_name}_!{chrom}_sample_metadata.tsv > modern_arc_groups.tsv
   bcftools +split -G modern_arc_groups.tsv -Oz -o . !{params.run_name}_!{chrom}.vcf.gz 2> bcftools_split_!{params.run_name}_!{chrom}.stderr > bcftools_split_!{params.run_name}_!{chrom}.stdout
   tabix -f -p vcf !{params.run_name}_!{chrom}_archaic.vcf.gz
   tabix -f -p vcf !{params.run_name}_!{chrom}_modern.vcf.gz
   '''
}

process vcfsubset {
   tag "${pop} chr${chrom}"

   cpus params.vcfsubset_cpus
   memory { params.vcfsubset_mem.plus(task.attempt.minus(1).multiply(params.vcfsubset_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.vcfsubset_retry_timeout : params.vcfsubset_timeout }
   queue { task.attempt >= 2 ? params.vcfsubset_retry_queue : params.vcfsubset_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*_samples.tsv', saveAs: { 'chr'+chrom+'_'+it }

   input:
   tuple val(pop), val(chrom), path(input_vcf), path(input_tbi), path(metadata)

   output:
   tuple path("bcftools_view_selectpops_Sprime_${pop}_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_Sprime_${pop}_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_Sprime_${pop}_chr${chrom}.stdout"), emit: logs
   tuple val(pop), path("${pop}_chr${chrom}_nomissinggenos.vcf.gz"), path("${pop}_chr${chrom}_nomissinggenos.vcf.gz.tbi"), emit: vcf
   tuple val(pop), path("${pop}_outgroup_samples.tsv"), emit: outgroup
   tuple val(pop), path("${pop}_samples.tsv"), emit: samples

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_outgroup_colname}" -v "select=!{params.sprime_outgroup}" !{metadata} > !{pop}_outgroup_samples.tsv
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{pop}" !{metadata} > !{pop}_samples.tsv
   cat !{pop}_samples.tsv !{pop}_outgroup_samples.tsv > !{pop}_plusOutgroup_samples.tsv
   #We have to separate the sample selection from the genotype missingness filter
   # because of a bug in bcftools view when -S and -g are used in tandem:
   #Way too many sites get filtered if you use them together, possibly more than
   # if -g was used alone.
   bcftools view -S !{pop}_plusOutgroup_samples.tsv -Ou !{input_vcf} 2> bcftools_view_selectpops_Sprime_!{pop}_chr!{chrom}.stderr | \
   bcftools view -g ^miss -Oz -o !{pop}_chr!{chrom}_nomissinggenos.vcf.gz 2> bcftools_view_nomissinggenos_Sprime_!{pop}_chr!{chrom}.stderr > bcftools_view_nomissinggenos_Sprime_!{pop}_chr!{chrom}.stdout
   tabix -f -p vcf !{pop}_chr!{chrom}_nomissinggenos.vcf.gz
   '''
}

process concatvcf {
   tag "${pop}"

   cpus params.concatvcf_cpus
   memory { params.concatvcf_mem.plus(task.attempt.minus(1).multiply(params.concatvcf_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.concatvcf_retry_timeout : params.concatvcf_timeout }
   queue { task.attempt >= 2 ? params.concatvcf_retry_queue : params.concatvcf_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), path("*"), path("*")
   val(autosome_num_list)
   val(num_autosomes)

   output:
   tuple path("bcftools_concat_Sprime_${pop}.stderr"), path("bcftools_concat_Sprime_${pop}.stdout"), emit: logs
   tuple val(pop), path("${pop}_nomissinggenos.vcf.gz"), path("${pop}_nomissinggenos.vcf.gz.tbi"), emit: vcf

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   if [[ "!{num_autosomes}" -gt "1" ]]; then
      bcftools concat --threads !{task.cpus} -Oz -o !{pop}_nomissinggenos.vcf.gz !{pop}_chr{!{autosome_num_list}}_nomissinggenos.vcf.gz 2> bcftools_concat_Sprime_!{pop}.stderr > bcftools_concat_Sprime_!{pop}.stdout
   else
      ln -s !{pop}_chr!{autosome_num_list}_nomissinggenos.vcf.gz !{pop}_nomissinggenos.vcf.gz
      echo "Only one chromosome provided, skipping bcftools concat for !{pop}" > bcftools_concat_Sprime_!{pop}.stderr
      echo "" > bcftools_concat_Sprime_!{pop}.stdout
   fi
   tabix -f -p vcf !{pop}_nomissinggenos.vcf.gz
   '''
}

process sprime {
   tag "${pop} chr${chrom}"

   cpus params.sprime_cpus
   memory { params.sprime_mem.plus(1).plus(task.attempt.minus(1).multiply(params.sprime_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.sprime_retry_timeout : params.sprime_timeout }
   queue { task.attempt >= 2 ? params.sprime_retry_queue : params.sprime_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.log'
   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.score'

   input:
   tuple val(pop), path("${pop}_nomissinggenos.vcf.gz"), path("${pop}_nomissinggenos.vcf.gz.tbi"), path("${pop}_outgroup_samples.tsv"), val(chrom), path(recmap)
   val(autosome_map)

   output:
   tuple path("Sprime_${pop}_chr${chrom}.stderr"), path("Sprime_${pop}_chr${chrom}.stdout"), path("${pop}_chr${chrom}_Sprime.log"), emit: logs
   tuple val(pop), path("${pop}_chr${chrom}_Sprime.score"), emit: score

   shell:
   sprime_retry_mem = params.sprime_mem.plus(task.attempt.minus(1).multiply(params.sprime_mem_increment))
   chrid = autosome_map[chrom]
   '''
   module load !{params.mod_Sprime}
   java -Xmx!{sprime_retry_mem}g -jar ${SPRIME} gt=!{pop}_nomissinggenos.vcf.gz outgroup=!{pop}_outgroup_samples.tsv map=!{recmap} chrom=!{chrid} out=!{pop}_chr!{chrom}_Sprime 2> Sprime_!{pop}_chr!{chrom}.stderr > Sprime_!{pop}_chr!{chrom}.stdout
   '''
}

process matchrate {
   tag "${pop}"

   cpus params.matchrate_cpus
   memory { params.matchrate_mem.plus(task.attempt.minus(1).multiply(params.matchrate_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.matchrate_retry_timeout : params.matchrate_timeout }
   queue { task.attempt >= 2 ? params.matchrate_retry_queue : params.matchrate_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.tsv.gz'

   input:
//This path("*") corresponds to the Sprime .score files:
//The following path("*") corresponds to the archaic VCFs, tabix indices, and
// the metadata file including archaic samples:
   tuple val(pop), path(scores), path(vcfs), path(tbis), path(metadatas)
   val(autosome_list)
   val(has_chr_prefix)

   output:
   tuple val(pop), path("${pop}_autosomes_Sprime_matches.tsv.gz"), emit: matches
   tuple val(pop), path("${pop}_autosomes_Sprime_match_rates.tsv.gz"), emit: matchrates

   shell:
   '''
   module load !{params.mod_bcftools}
   #Prep the archaic metadata file from the per-chrom sample metadata files:
   #Columns of archaic metadata file must be named Sample and Region, and
   # the Regions should be Neandertal (and Denisovan, if included) for
   # the match rate scripts to work, rather than e.g. NeaA and DenA.
   !{projectDir}/HumanPopGenScripts/Simulations/prep_arc_metadata.awk -v "simgroups=!{params.sim_arc_groups}" -v "groups=!{params.archaic_groups}" *_sample_metadata.tsv > archaic_metadata.tsv
   #Identify Sprime alleles that match or mismatch the archaics:
   #The extra comma is not a typo, it's a way to handle only one autosome with brace expansion.
   for chr in {!{autosome_list},};
      do
      if [[ "!{has_chr_prefix}" == "true" ]]; then
         chrnum=${chr#chr};
      else
         chrnum=${chr};
      fi;
      #Establish the set of Sprime sites for this population and chromosome:
      fgrep -v "CHROM" !{pop}_chr${chrnum}_Sprime.score | \
         cut -f1,2 > !{pop}_chr${chrnum}_Sprime_sites.tsv
      !{projectDir}/HumanPopGenScripts/Sprime/archaicMatchSprime.awk -v "spop=!{pop}" -v "groups=!{params.archaic_groups}" -v "arcs=!{params.archaic_indivs}" -v "short=!{params.archaic_shortnames}" archaic_metadata.tsv \
         <(bcftools view -R !{pop}_chr${chrnum}_Sprime_sites.tsv -T !{pop}_chr${chrnum}_Sprime_sites.tsv !{params.run_name}_${chr}_archaic.vcf.gz) \
         !{pop}_chr${chrnum}_Sprime.score
   done | \
      awk 'BEGIN{FS="\t";OFS=FS;}/^CHROM/&&NR==1{print;}!/^CHROM/{print;}' | \
      gzip -9 > !{pop}_autosomes_Sprime_matches.tsv.gz
   #Calculate the match rates from this file:
   gzip -dc !{pop}_autosomes_Sprime_matches.tsv.gz | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeArchaicMatchRate.awk -v "species=!{params.archaic_groups}" -v "arcs=!{params.archaic_shortnames}" | \
      gzip -9 > !{pop}_autosomes_Sprime_match_rates.tsv.gz
   '''
}

process projection {
   tag "${pop}"

   cpus params.projection_cpus
   memory { params.projection_mem.plus(task.attempt.minus(1).multiply(params.projection_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.projection_retry_timeout : params.projection_timeout }
   queue { task.attempt >= 2 ? params.projection_retry_queue : params.projection_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.bed'

   input:
//These path("*") correspond to Sprime .score files, modern VCFs, tbis,
// and outgroup sample lists, respectively:
   tuple val(pop), path(scores), path(vcfs), path(tbis), path(outgroupsamples)
   val(autosome_list)
   val(has_chr_prefix)

   output:
   tuple val(pop), path("${pop}_Sprime_phased_tracts_perSample_maxgap${params.tract_max_gap}.bed"), emit: projections

   shell:
   '''
   module load !{params.mod_bcftools}
   #Project Sprime alleles onto each individual to identify tracts from genotypes:
   #The extra comma is not a typo, it's a way to handle only one autosome with brace expansion.
   for chr in {!{autosome_list},};
      do
      if [[ "!{has_chr_prefix}" == "true" ]]; then
         chrnum=${chr#chr};
      else
         chrnum=${chr};
      fi;
      bcftools query -f '%CHROM:%POS[\t%GT]\n' -r ${chr} -H -S ^!{pop}_outgroup_samples.tsv !{pop}_chr${chrnum}_nomissinggenos.vcf.gz | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimePerSampleTracts.awk -v "spop=!{pop}" -v "allout=1" -v "phased=1" !{pop}_chr${chrnum}_Sprime.score - | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBED.awk -v "phased=1" -v "max_gap=!{params.tract_max_gap}" | \
         sort -k1,1V -k2,2n -k3,3n
   done > !{pop}_Sprime_phased_tracts_perSample_maxgap!{params.tract_max_gap}.bed
   '''
}

process tagsites {
   tag "${pop} ${origin}"

   cpus params.tagsites_cpus
   memory { params.tagsites_mem.plus(task.attempt.minus(1).multiply(params.tagsites_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.tagsites_retry_timeout : params.tagsites_timeout }
   queue { task.attempt >= 2 ? params.tagsites_retry_queue : params.tagsites_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(pop), path(matches), path(matchrates), val(origin)

   output:
   tuple val(pop), val(origin), path("${pop}_${params.origin_criteria}_${origin}Match_sites.tsv"), emit: sites

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
   memory { params.calcld_mem.plus(task.attempt.minus(1).multiply(params.calcld_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.calcld_retry_timeout : params.calcld_timeout }
   queue { task.attempt >= 2 ? params.calcld_retry_queue : params.calcld_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), val(origin), path(sites), val(chrom), path(vcf), path(tbi), path(metadata)
   val(autosome_map)

   output:
   tuple path("bcftools_view_${pop}_chr${chrom}_${origin}.stderr"), path("vcftools_${pop}_chr${chrom}_${params.r2_type}ld_r${params.min_r2}_${origin}.stderr"), path("vcftools_${pop}_chr${chrom}_${params.r2_type}ld_r${params.min_r2}_${origin}.stdout"), emit: logs
   tuple val(pop), path("${pop}_chr${chrom}_Sprime_${origin}_${params.r2_type}ld_r${params.min_r2}_core_haplotypes.tsv.gz"), emit: corehaps

   shell:
   chrid = autosome_map[chrom]
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_vcftools}
   module load !{params.mod_R}
   #Identify the set of samples from the superpopulation containing the
   # target population, since we want LD calculated from a large sample:
   superpop=$(!{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.r2_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{pop}" !{metadata} | uniq)
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.r2_colname}" -v "select=${superpop}" !{metadata} > !{pop}_superpop_samples.tsv
   #Pull out the archaic-matching S' sites and calculate pairwise r^2:
   awk -v "chrom=!{chrid}" '\$1==chrom' !{sites} > !{origin}_sites.tsv
   if [[ -s "!{pop}_superpop_samples.tsv" ]] && [[ -s "!{origin}_sites.tsv" ]]; then
      bcftools view -R !{origin}_sites.tsv -T !{origin}_sites.tsv -S !{pop}_superpop_samples.tsv -a -Ov !{vcf} 2> bcftools_view_!{pop}_chr!{chrom}_!{origin}.stderr | \
         vcftools --vcf - --out !{pop}_chr!{chrom}_!{origin}_LD_trimalt --!{params.r2_type}-r2 --min-r2 !{params.min_r2} 2> vcftools_!{pop}_chr!{chrom}_!{params.r2_type}ld_r!{params.min_r2}_!{origin}.stderr > vcftools_!{pop}_chr!{chrom}_!{params.r2_type}ld_r!{params.min_r2}_!{origin}.stdout
      #Filter and reformat the r^2 results to only between pairs of input sites:
      !{projectDir}/HumanPopGenScripts/Sprime/extract_arcmatch_tag_SNPs_Rsquared.awk -v "source=vcftools" !{origin}_sites.tsv !{pop}_chr!{chrom}_!{origin}_LD_trimalt.!{params.r2_type}.ld > !{pop}_chr!{chrom}_!{origin}_!{params.r2_type}ld_r!{params.min_r2}_rsquared.tsv
      #Construct a graph from these r^2 edges and identify core haplotypes as
      # connected components of the graph, outputting a labeled list of sites:
      !{projectDir}/HumanPopGenScripts/Sprime/LD_connected_components.R !{pop}_chr!{chrom}_!{origin}_!{params.r2_type}ld_r!{params.min_r2}_rsquared.tsv !{pop}_chr!{chrom}_Sprime_!{origin}_!{params.r2_type}ld_r!{params.min_r2}_core_haplotypes.tsv.gz !{chrid}
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

process catcorehaps {
   tag "${pop}"

   cpus params.catcorehaps_cpus
   memory { params.catcorehaps_mem.plus(task.attempt.minus(1).multiply(params.catcorehaps_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.catcorehaps_retry_timeout : params.catcorehaps_timeout }
   queue { task.attempt >= 2 ? params.catcorehaps_retry_queue : params.catcorehaps_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), path(corehaps)

   output:
   tuple val(pop), path("${pop}_Sprime_${params.r2_type}ld_r${params.min_r2}_core_haplotypes_annotated.tsv.gz"), emit: corehaps

   shell:
   '''
   #Stage the core haplotype files:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[0]} ${a[1]};
   done < !{corehaps}
   #Concatenate the core haplotypes across both chromosomes and tract origins:
   while read a;
      do
      fn=$(basename ${a});
      prefix=${fn%_!{params.r2_type}ld_r!{params.min_r2}_core_haplotypes.tsv.gz};
      origin=${prefix#!{pop}_chr*_Sprime_};
      !{projectDir}/HumanPopGenScripts/Sprime/cat_core_haplotypes.awk -v "pop=!{pop}" -v "origin=${origin}" <(gzip -dc ${a});
   done < <(cut -f2 !{corehaps} | sort -t"_" -k1,1 -k2,2V -k4,4) | \
      gzip -9 > !{pop}_Sprime_!{params.r2_type}ld_r!{params.min_r2}_core_haplotypes_annotated.tsv.gz
   '''
}

process core_score {
   tag "${pop}"

   cpus params.corescore_cpus
   memory { params.corescore_mem.plus(task.attempt.minus(1).multiply(params.corescore_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.corescore_retry_timeout : params.corescore_timeout }
   queue { task.attempt >= 2 ? params.corescore_retry_queue : params.corescore_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/corehaps", mode: 'copy', pattern: '*.tsv'

   input:
   tuple val(pop), path(corehaps), path(matches)

   output:
   tuple val(pop), path("${pop}_Sprime_corehaps_matches.tsv"), emit: matches

   shell:
   '''
   !{projectDir}/HumanPopGenScripts/Sprime/applyCoreHaplotypesToScore.awk <(gzip -dc !{corehaps}) <(gzip -dc !{matches}) > !{pop}_Sprime_corehaps_matches.tsv
   '''
}

process corehapfreqs {
   tag "${pop} chr${chrom}"

   cpus params.coretf_cpus
   memory { params.coretf_mem.plus(task.attempt.minus(1).multiply(params.coretf_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.coretf_retry_timeout : params.coretf_timeout }
   queue { task.attempt >= 2 ? params.coretf_retry_queue : params.coretf_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   tuple val(pop), path(corehapscores), val(chrom), path(vcf), path(tbi), path(metadata)
   val(autosome_map)

   output:
   tuple val(pop), val(chrom), path("${pop}_chr${chrom}_Sprime_${pop}_corehap_freqs.tsv.gz"), emit: freqs

   shell:
   chrid = autosome_map[chrom]
   '''
   module load !{params.mod_bcftools}
   #Identify the samples for this population:
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{pop}" !{metadata} > !{pop}_pop_samples.tsv
   #Calculate the median Sprime allele frequency for each tract:
   bcftools query -f '%CHROM:%POS[\t%GT]\n' -r !{chrid} -H -S !{pop}_pop_samples.tsv !{vcf} | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeArchaicAF.awk -v "spop=!{pop}" -v "pop=!{pop}" -v "all=1" !{corehapscores} - | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractMedianAF.awk | \
      gzip -9 > !{pop}_chr!{chrom}_Sprime_!{pop}_corehap_freqs.tsv.gz
   '''
}

process tract_map {
   cpus params.tractmap_cpus
   memory { params.tractmap_mem.plus(task.attempt.minus(1).multiply(params.tractmap_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.tractmap_retry_timeout : params.tractmap_timeout }
   queue { task.attempt >= 2 ? params.tractmap_retry_queue : params.tractmap_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/corehaps", mode: 'copy', pattern: '*.tsv.gz'

   input:
   tuple path(corescores), path(targetpops)

   output:
   path "${params.run_name}_Sprime_corehap_TractID_map.tsv.gz", emit: tractmap

   shell:
   '''
   #Symlink the many input files to avoid SLURM SBATCH script size limits:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[2]} ${a[1]};
   done < !{corescores}
   #Generate a map of core haplotype IDs to tract IDs from the adjusted
   # matches/.score files:
   header=1
   while IFS=$'\t' read p;
      do
      !{projectDir}/HumanPopGenScripts/Sprime/corehapTractIDmap.awk -v "header=${header}" ${p}_Sprime_corehaps_matches.tsv
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_corehap_TractID_map.tsv.gz
   unset header
   '''
}

process pbs {
   tag "${target} ${chrom}"

   cpus params.pbs_cpus
   memory { params.pbs_mem.plus(task.attempt.minus(1).multiply(params.pbs_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.pbs_retry_timeout : params.pbs_timeout }
   queue { task.attempt >= 2 ? params.pbs_retry_queue : params.pbs_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(vcf), path(tbi), path(metadata), val(target), val(ref), val(og)

   output:
   tuple path("pbs_cli_${params.run_name}_${target}_chr${chrom}_w${params.window_size}_step${params.window_step}.stderr"), path("pbs_cli_${params.run_name}_${target}_chr${chrom}_w${params.window_size}_step${params.window_step}.stdout"), emit: logs
   tuple val(target), val(chrom), path("${params.run_name}_${target}_chr${chrom}_pbs_w${params.window_size}_step${params.window_step}.tsv.gz"), emit: scores

   shell:
   '''
   module load !{params.mod_miniconda}
   conda activate adintsims
   !{projectDir}/pbs_cli.py -i !{vcf} -m !{metadata} -s SampleID -p Population --superpop_col Region -t !{target} -r !{ref} -e !{og} -w !{params.window_size} --window_step !{params.window_step} !{params.clip_pbs} -o !{params.run_name}_!{target}_chr!{chrom}_pbs_w!{params.window_size}_step!{params.window_step}.tsv.gz 2> pbs_cli_!{params.run_name}_!{target}_chr!{chrom}_w!{params.window_size}_step!{params.window_step}.stderr > pbs_cli_!{params.run_name}_!{target}_chr!{chrom}_w!{params.window_size}_step!{params.window_step}.stdout
   '''
}

process xpehh {
   tag "${target} ${chrom}"

   cpus params.xpehh_cpus
   memory { params.xpehh_mem.plus(task.attempt.minus(1).multiply(params.xpehh_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.xpehh_retry_timeout : params.xpehh_timeout }
   queue { task.attempt >= 2 ? params.xpehh_retry_queue : params.xpehh_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(vcf), path(tbi), path(metadata), path(gmap), val(target), val(ref)

   output:
   tuple path("xpehh_cli_${params.run_name}_${target}_chr${chrom}.stderr"), path("xpehh_cli_${params.run_name}_${target}_chr${chrom}.stdout"), emit: logs
   tuple val(target), val(chrom), path("${params.run_name}_${target}_chr${chrom}_xpehh.tsv.gz"), emit: scores

   shell:
   '''
   module load !{params.mod_miniconda}
   conda activate adintsims
   OMP_NUM_THREADS=!{task.cpus} !{projectDir}/xpehh_cli.py -i !{vcf} -m !{metadata} -s SampleID -p Population --superpop_col Region -g !{gmap} --pos_col_index !{params.gmap_pos_col} --cm_col_index !{params.gmap_cm_col} -t !{target} -r !{ref} -o !{params.run_name}_!{target}_chr!{chrom}_xpehh.tsv.gz 2> xpehh_cli_!{params.run_name}_!{target}_chr!{chrom}.stderr > xpehh_cli_!{params.run_name}_!{target}_chr!{chrom}.stdout
   '''
}

process cat_outs {
   cpus params.catouts_cpus
   memory { params.catouts_mem.plus(task.attempt.minus(1).multiply(params.catouts_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.catouts_retry_timeout : params.catouts_timeout }
   queue { task.attempt >= 2 ? params.catouts_retry_queue : params.catouts_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.t{sv,csv}.gz'

   input:
   tuple path(matchrates), path(projections), path(tractfreqs), path(pbs), path(xpehh), path(targetpops), path(metadata)
   val(autosome_num_list)

   output:
   tuple path("${params.run_name}_perPop_Sprime_autosomal_match_rates.tsv.gz"), path("${params.run_name}_Sprime_perChrom_perIndiv_perPop_tract_lengths.tsv.gz"), path("${params.run_name}_Sprime_targetpop_corehap_freqs.tsv.gz"), path("${params.run_name}_perPop_pbs_w${params.window_size}_step${params.window_step}.tsv.gz"), path("${params.run_name}_perPop_xpehh.tsv.gz"), emit: outputs

   shell:
   '''
   #Set the array of chromosomes to iterate over:
   IFS="," read -a chroms <<< "!{autosome_num_list}"
   #Symlink the many input files to avoid SLURM SBATCH script size limits:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[2]} ${a[1]};
   done < <(cat !{matchrates} !{projections} !{tractfreqs} !{pbs} !{xpehh})
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
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBEDtoLengths.awk -v "header=${header}" -v "pop=${p}" -v "phased=1" -v "addpop=0" - ${p}_Sprime_phased_tracts_perSample_maxgap!{params.tract_max_gap}.bed
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_perChrom_perIndiv_perPop_tract_lengths.tsv.gz
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
   #PBS:
   header=1
   while IFS=$'\t' read p;
      do
      for chrom in ${chroms[@]};
         do
         gzip -dc !{params.run_name}_${p}_chr${chrom}_pbs_w!{params.window_size}_step!{params.window_step}.tsv.gz | \
            awk -v "header=${header}" 'BEGIN{FS="\t";OFS=FS;}NR==1&&header{print;}NR>1{print;}'
         header=0;
      done;
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_perPop_pbs_w!{params.window_size}_step!{params.window_step}.tsv.gz
   unset header
   #XP-EHH:
   header=1
   while IFS=$'\t' read p;
      do
      for chrom in ${chroms[@]};
         do
         gzip -dc !{params.run_name}_${p}_chr${chrom}_xpehh.tsv.gz | \
            awk -v "header=${header}" 'BEGIN{FS="\t";OFS=FS;}NR==1&&header{print;}NR>1{print;}'
         header=0;
      done;
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_perPop_xpehh.tsv.gz
   unset header
   '''
}

workflow {
   //Detection mechanism for "chr" prefixes of autosomes:
   has_chr_prefix = params.autosomes.count('chr') > 0
   autosome_list = params.autosomes
   autosome_num_list = has_chr_prefix ? params.autosomes.replaceAll('chr', '') : autosome_list
   Map autosome_map = [params.autosomes.replaceAll('chr', '').tokenize(','), params.autosomes.tokenize(',')]
      .transpose()
      .collectEntries()

   num_autosomes = params.autosomes.tokenize(',').size()

   //Create channel of the autosomes for simulations:
   Channel
      .fromList(autosome_num_list.tokenize(','))
      .tap { autosomes }
      .subscribe { println "Added ${it} to autosomes channel" }

   //VCFs get generated by the simulations, so no need to read them in

   //Set up the channel of recombination rate map files for Sprime:
   Channel
      .fromPath(params.recmap_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find recombination rate map files matching glob: ${params.recmap_glob}" }
      .map { a -> [ (a =~ params.recmap_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from recombination rate map: ${params.recmap_glob} =~ ${params.recmap_regex}" }
      .tap { recmaps }
      .subscribe { println "Added chr${it[0]} recombination rate map (${it[1]}) to recmaps channel" }

   //Retain only the autosomes for the recombination rate maps:
   recmaps_autosomes = recmaps
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
      .tap { recmaps_autosomes_for_xpehh }

   //Set up the channel of target populations for Sprime:
   Channel
      .fromPath(params.target_pops_file, checkIfExists: true)
      .ifEmpty { error "Unable to find target populations file for Sprime, PBS, and XP-EHH: ${params.target_pops_file}" }
      .splitCsv(sep: "\t")
      .tap { pbs_pops }
      .map { target,ref,og -> [ target, ref ] }
      .tap { xpehh_pops }
      .map { target,ref -> target }
      .tap { sprime_pops }
      .tap { sprime_target_pops }
      .tap { tractmap_target_pops }
      .subscribe { println "Added ${it} to *_pops channels" }

   //Set up the channel of tract origins:
   Channel
      .fromList(params.tract_origin_list.tokenize(','))
      .tap { tract_origins }
      .subscribe { println "Added ${it} to tract_origins channel" }

   //Run the simulations for each chromosome, convert to VCF, and split
   // into archaics-only and moderns-only VCFs:
   n_autosomes = Channel.value(num_autosomes)
   run_sim(autosomes, n_autosomes, autosome_map)

   //Subset the input VCFs for use with Sprime (i.e. retain target and outgroup pops and no missing genotypes):
   vcfsubset_inputs = sprime_pops
      .combine(run_sim.out.modern_vcf
                  .join(run_sim.out.metadata, by: 0, failOnDuplicate: true, failOnMismatch: true))
   vcfsubset(vcfsubset_inputs)
//vcfsubset.out.logs into sprime_vcf_subset_logs
//vcfsubset.out.vcf into Sprime_perchrom_vcfs,Sprime_perchrom_vcfs_project,Sprime_perchrom_vcfs_tractfreqs
//vcfsubset.out.outgroup into Sprime_outgroup,Sprime_outgroup_project
//vcfsubset.out.samples into Sprime_pop_samples

   //Prepare VCFs for concatvcf and projection steps:
   vcfs_for_concatvcf = vcfsubset.out.vcf
      .groupTuple(by: 0, size: num_autosomes)
      .tap { vcfs_for_projection }

   //Prepare single outgroup files for Sprime and projection steps:
   outgroups_for_sprime = vcfsubset.out.outgroup
      .groupTuple(by: 0, size: num_autosomes)
      .map({ a -> [ a[0], a[1].take(1) ] })
      .tap { outgroups_for_projection }

   //Concatenate the per-chromosome subsetted VCFs to serve as Sprime input:
   concatvcf(vcfs_for_concatvcf, autosome_num_list, num_autosomes)
//concatvcf.out.logs into concat_input_vcf_logs
//concatvcf.out.vcf into Sprime_vcfs

   //Run Sprime per chromosome but with concatenated VCF:
   sprime_inputs = concatvcf.out.vcf
      .join(outgroups_for_sprime, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .combine(recmaps_autosomes)
   sprime(sprime_inputs, autosome_map)
//sprime.out.log into Sprime_logs
//sprime.out.score into Sprime_scores_matchrate,Sprime_scores_project,Sprime_scores_tractfreqs

   //Calculate the match rates of each tract to the high-depth archaics:
   matchrate_inputs = sprime.out.score
      .groupTuple(by: 0, size: num_autosomes)
      .combine(run_sim.out.archaic_vcf
         .join(run_sim.out.metadata, by: 0, failOnDuplicate: true, failOnMismatch: true)
         .map({ chrom,vcf,tbi,metadata -> [ "arc", vcf, tbi, metadata ] })
         .groupTuple(by: 0, size: num_autosomes)
         .map({ arc,vcfs,tbis,metadatas -> [ vcfs, tbis, metadatas ] }))
   matchrate(matchrate_inputs, autosome_list, has_chr_prefix)
//matchrate.out.matches into Sprime_matches
//matchrate.out.matchrates into Sprime_matchrates

   //Project Sprime population-level tracts onto haplotypes:
   projection_inputs = sprime.out.score
      .groupTuple(by: 0, size: num_autosomes)
      .join(vcfsubset.out.vcf.groupTuple(by: 0, size: num_autosomes), by: 0, failOnDuplicate: true, failOnMismatch: true)
      .join(outgroups_for_projection, by: 0, failOnDuplicate: true, failOnMismatch: true)
   projection(projection_inputs, autosome_list, has_chr_prefix)
//projection.out.projections into Sprime_project_BEDs

   //Identify Sprime "tag" sites:
   tagsites_inputs = matchrate.out.matches
      .join(matchrate.out.matchrates, by: 0, failOnMismatch: true, failOnDuplicate: true)
      .combine(tract_origins)
   tagsites(tagsites_inputs)
//tuple val(pop), path(matches), path(matchrates), val(origin)
//tagsites.out.sites into tagsites_output

   //Calculate LD amongst tag sites and identify core haplotypes:
   calc_ld_inputs = tagsites.out.sites
      .combine(run_sim.out.modern_vcf.join(run_sim.out.metadata, by: 0, failOnMismatch: true, failOnDuplicate: true))
   calc_ld(calc_ld_inputs, autosome_map)
//tuple val(pop), val(origin), path(sites), val(chrom), path(vcf), path(tbi), path(metadata)
//calc_ld.out.logs into calc_ld_logs
//calc_ld.out.corehaps into calc_ld_output

   //Concatenate core haplotypes across chroms and tract origins:
   catcorehaps_inputs = calc_ld.out.corehaps
      .collectFile() { a -> [ "${a[0]}.txt", a[1].toAbsolutePath().toString()+'\t'+a[1].getName()+'\n' ] }
      .map({ a -> [ a.getSimpleName(), a ] })
   catcorehaps(catcorehaps_inputs)
//tuple val(pop), path(corehaps)
//catcorehaps.out.corehaps into cat_corehaps_output

   //Annotate matches file with core haplotype IDs:
   core_score_inputs = catcorehaps.out.corehaps
      .join(matchrate.out.matches, by: 0, failOnDuplicate: true, failOnMismatch: true)
   core_score(core_score_inputs)
//tuple val(pop), path(corehaps), path(matches)
//core_score.out.matches into corehaps_score_for_project,corehaps_score_for_tractfreqs,corehaps_score_for_tractmap

   //Estimate the frequency of each core haplotype:
   corehapfreqs_inputs = core_score.out.matches
      .combine(run_sim.out.modern_vcf.join(run_sim.out.metadata, by: 0, failOnDuplicate: true, failOnMismatch: true))
   corehapfreqs(corehapfreqs_inputs, autosome_map)
//tuple val(pop), path(corehapscores), val(chrom), path(vcf), path(tbi), path(metadata)
//corehapfreqs.out.freqs into corehaps_tract_freqs_for_genes,corehaps_tract_freqs_for_cat

   //Generate a map from core haplotype IDs to Sprime tract IDs:
   tract_map_inputs = core_score.out.matches
      .collectFile() { pop,matches -> [ "score_paths.tsv", matches.getSimpleName()+'\t'+matches.getName()+'\t'+matches+'\n' ] }
      .combine(tractmap_target_pops.collectFile(name: 'Sprime_target_populations.txt', newLine: true, sort: true))
   tract_map(tract_map_inputs)
//tuple path(corescores), path(targetpops)
//tract_map.out.tractmap into tract_map_output

   //Run PBS:
   pbs_inputs = run_sim.out.modern_vcf
      .join(run_sim.out.metadata, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .combine(pbs_pops)
   pbs(pbs_inputs)

   //Run XP-EHH:
   xpehh_inputs = run_sim.out.modern_vcf
      .join(run_sim.out.metadata, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .join(recmaps_autosomes_for_xpehh, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .combine(xpehh_pops)
   xpehh(xpehh_inputs)

   //Concatenate outputs for simplicity of analysis:
   cat_outs_inputs = matchrate.out.matchrates
      .collectFile() { pop,matchrates -> [ "matchrate_paths.tsv", matchrates.getSimpleName()+'\t'+matchrates.getName()+'\t'+matchrates+'\n' ] }
      .combine(projection.out.projections.collectFile() { pop,projections -> [ "projection_paths.tsv", projections.getSimpleName()+'\t'+projections.getName()+'\t'+projections+'\n' ] })
      .combine(corehapfreqs.out.freqs.collectFile() { pop,chr,freqs -> [ "corehapfreq_paths.tsv", freqs.getSimpleName()+'\t'+freqs.getName()+'\t'+freqs+'\n' ] })
      .combine(pbs.out.scores.collectFile() { pop,chr,pbs -> [ "pbs_paths.tsv", pbs.getSimpleName()+'\t'+pbs.getName()+'\t'+pbs+'\n' ] })
      .combine(xpehh.out.scores.collectFile() { pop,chr,xpehh -> [ "xpehh_paths.tsv", xpehh.getSimpleName()+'\t'+xpehh.getName()+'\t'+xpehh+'\n' ] })
      .combine(sprime_target_pops.collectFile(name: 'Sprime_target_populations.txt', newLine: true, sort: true))
      .combine(run_sim.out.metadata.map({ chr,meta -> meta }).first())
   cat_outs(cat_outs_inputs, autosome_num_list)
//cat_outs.out.outputs into catouts_outputs
}

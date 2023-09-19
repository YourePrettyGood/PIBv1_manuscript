#!/usr/bin/env nextflow
/* Pipeline to generate haplotype plots and NEXUS files for   *
 * haplotype networks from genomic regions containing archaic *
 * introgressed tracts.                                       *
 * Core steps:                                                *
 *  Extract sites with MAF >= threshold in modern samples ->  *
 *  Extract the same sites from archaic samples ->            *
 *  Merge VCFs -> Add ancestral allele INFO tags to SNPs ->   *
 *  |>Convert haplotypes to TSV for R and NEXUS for networks  *
 *  |>Extract CADD scores for SNPs                            *
 *  ->Run VEP and extract the worst consequences              *
 * Extra steps:                                               *
 *  Extract MPRA emVars in the region                         *
 *  Extract Sprime tract and core haplotype sites             *
 *  Convert rsids of published variants to asm coordinates    *
 *  Construct a matrix of annotated site presence             */

nextflow.enable.dsl=2

//Default paths, globs, and regexes:
//Modern sample VCFs (phased):
params.modern_glob = "${projectDir}/modern_VCFs/*.vcf.gz"
params.modern_regex = ~/_chr(\p{Alnum}+)_phased$/
//Archaic VCFs (phased):
params.archaic_glob = "${projectDir}/archaic_VCFs/*.vcf.gz"
params.archaic_regex = ~/_chr(\p{Alnum}+)_pop_phased$/
//Regions to plot, in UCSC-like coordinates (i.e. [chromosome]:[1-based start]-[1-based end]):
//This file should consist of seven columns separated by tabs:
// 1) Chromosome of the region
// 2) Short name (alphanumeric, no special characters) of the region (used for filenames)
// 3) Region coordinates in UCSC-like format
// 4) Sprime tract ID
// 5) Sprime core haplotype ID
// 6) Pilot Sprime tract ID
// 7) Path to file of paths to rsids of published variants (can be empty file, files of rsids should be named *_[locus short name].txt)
params.target_regions = "${projectDir}/${params.run_name}_target_regions.tsv"
//Minimum global MAF for filtering:
params.minmaf = "0.01"
//Samples to keep for plotting:
params.samples_to_keep = "${projectDir}/${params.run_name}_samples_to_plot.txt"
//Ancestral allele FASTAs:
params.aa_glob = "/gpfs/gibbs/pi/tucci/pfr8/human_ancestor/Ensembl_R75_GRCh37/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_*.fa"
params.aa_regex = ~/ancestor_(\p{Alnum}+)$/
//TSVs from WGSA of CADD 1.6 scores:
params.cadd_glob = "/gpfs/gibbs/pi/tucci/pfr8/WGSA/resources/CADDv1.6/hg19/whole_genome_SNVs.tsv.gz.chr*.gz.rankRawScore.gz"
params.cadd_regex = ~/[.]chr(\p{Alnum}+)[.]gz[.]rankRawScore[.]gz$/
//GFF3 annotation of the reference:
params.ref_annotation = "/gpfs/gibbs/pi/tucci/pfr8/refs/GENCODE/r38/gencode.v38lift37.annotation.gff3.gz"
//MPRA emVar BEDs:
params.k562_bed = "/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/MPRA/PIB92_pilot/results/K562/adaptive_variants_MPRA_K562_variant_ucsc_custom_track.bed"
params.jurkat_bed = "/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/MPRA/PIB92_pilot/results/Jurkat/adaptive_variants_MPRA_Jurkat_variant_ucsc_custom_track.bed"
//Final Sprime score files:
params.final_sprime_glob = "/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes/PIBv1_Sprime/Sprime/*_Sprime.score"
params.final_sprime_regex = ~/^(\p{Alnum}+)_chr(\p{Alnum}+)_Sprime$/
//Core haplotype matches files:
params.corehap_glob = "/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes/PIBv1_Sprime/corehaps/*_corehaps_matches.tsv"
params.corehap_regex = ~/^(\p{Alnum}+)_Sprime_corehaps_matches$/
//Pilot Sprime score files:
params.pilot_glob = "/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/FullDepthBatches/First92_analysis/hs37d5/joint_genotyping_with_SerenaData/Sprime/pop_results/PIB92_SD_hs37d5_SNP_VQSR99.5_INDEL_VQSR99.0_dbSNP138_VQSRpass_chr*_Sprime_defaults.score"
params.pilot_regex = ~/_chr(\p{Alnum}+)_(\p{Alnum}+)_Sprime_defaults$/
//dbSNP VCF:
params.dbsnp_vcf = "/gpfs/gibbs/pi/tucci/pfr8/Archaic/ACVD_merged_dbSNP154/final_VCFs/dbSNP_154_hs37d5.vcf.gz"
//VEP parameters:
params.vep_cache = "/gpfs/gibbs/pi/tucci/pfr8/vep_r107_cache/"
params.vep_species = "homo_sapiens"
params.vep_assembly = "GRCh37"

//Sample-Region map including archaics:
params.samplemap = "${projectDir}/${params.run_name}_sample_region_map.tsv"
//Column names for sample map:
//Sample ID column name:
params.samplecol = "SampleID"
//Region column name:
params.regioncol = "Region"

//Defaults for cpus, base memory, memory ramp, memory unit, base timeout, timeout ramp, and timeout unit
//Modern samples
params.modern_cpus = 1
params.modern_mem = 1
params.modern_mem_ramp = 8
params.modern_mem_unit = 'GB'
params.modern_timeout = 1
params.modern_timeout_ramp = 4
params.modern_timeout_unit = 'h'
//Archaic hominins
params.archaic_cpus = 1
params.archaic_mem = 1
params.archaic_mem_ramp = 8
params.archaic_mem_unit = 'GB'
params.archaic_timeout = 1
params.archaic_timeout_ramp = 4
params.archaic_timeout_unit = 'h'
//Merge
params.merge_cpus = 1
params.merge_mem = 1
params.merge_mem_ramp = 8
params.merge_mem_unit = 'GB'
params.merge_timeout = 1
params.merge_timeout_ramp = 4
params.merge_timeout_unit = 'h'
//Add ancestral alleles
params.addanc_cpus = 1
params.addanc_mem = 1
params.addanc_mem_ramp = 8
params.addanc_mem_unit = 'GB'
params.addanc_timeout = 1
params.addanc_timeout_ramp = 4
params.addanc_timeout_unit = 'h'
//Output haplotypes
params.hapsout_cpus = 1
params.hapsout_mem = 1
params.hapsout_mem_ramp = 8
params.hapsout_mem_unit = 'GB'
params.hapsout_timeout = 1
params.hapsout_timeout_ramp = 4
params.hapsout_timeout_unit = 'h'
//Extract CADD scores
params.cadd_cpus = 1
params.cadd_mem = 1
params.cadd_mem_ramp = 8
params.cadd_mem_unit = 'GB'
params.cadd_timeout = 1
params.cadd_timeout_ramp = 4
params.cadd_timeout_unit = 'h'
//Overlapping PCGs
params.pcgs_cpus = 1
params.pcgs_mem = 1
params.pcgs_mem_ramp = 8
params.pcgs_mem_unit = 'GB'
params.pcgs_timeout = 1
params.pcgs_timeout_ramp = 4
params.pcgs_timeout_unit = 'h'
//VEP
params.vep_cpus = 1
params.vep_mem = 4
params.vep_mem_ramp = 8
params.vep_mem_unit = 'GB'
params.vep_timeout = 1
params.vep_timeout_ramp = 4
params.vep_timeout_unit = 'h'
//Site-presence matrix
params.matrix_cpus = 1
params.matrix_mem = 1
params.matrix_mem_ramp = 8
params.matrix_mem_unit = 'GB'
params.matrix_timeout = 1
params.matrix_timeout_ramp = 4
params.matrix_timeout_unit = 'h'

process modern {
   tag "${regionname}"

   cpus params.modern_cpus
   memory { params.modern_mem.plus(task.attempt.minus(1).multiply(params.modern_mem_ramp))+' '+params.modern_mem_unit }
   time { params.modern_timeout.plus(task.attempt.minus(1).multiply(params.modern_timeout_ramp))+params.modern_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/VCFs", mode: 'copy', pattern: '*.vcf.g{z,z.tbi}'
   publishDir path: "${params.output_dir}/annotations", mode: 'copy', pattern: '*.tsv'

   input:
   tuple val(regionname), val(region), path(modernvcf), path(moderntbi)
   path samples

   output:
   tuple val(regionname), path("${outprefix}_MAFge${params.minmaf}.vcf.gz"), path("${outprefix}_MAFge${params.minmaf}.vcf.gz.tbi"), emit: mafvcf
   tuple val(regionname), path("${params.run_name}_${regionname}_modernMAFfiltered_sites.tsv"), emit: mafsites

   shell:
   outprefix="${params.run_name}_moderns_${regionname}_regionspan"
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   modernregion="!{outprefix}.vcf.gz"
   echo "Extracting sites in !{region} for samples (!{samples}) from !{modernvcf} into ${modernregion}"
   bcftools view -r !{region} -S !{samples} -a -Oz -o ${modernregion} !{modernvcf}
   tabix -f ${modernregion}
   modernregionminmaf="!{outprefix}_MAFge!{params.minmaf}.vcf.gz"
   minmafsites="!{params.run_name}_!{regionname}_modernMAFfiltered_sites.tsv"
   echo "Keeping only sites with nonmajor AF >= !{params.minmaf} from ${modernregion} into ${modernregionminmaf}"
   bcftools view -q !{params.minmaf}:nonmajor -Oz -o ${modernregionminmaf} ${modernregion}
   tabix -f ${modernregionminmaf}
   bcftools query -f '%CHROM\\t%POS\\n' ${modernregionminmaf} > ${minmafsites}
   '''
}

process archaic {
   tag "${regionname}"

   cpus params.archaic_cpus
   memory { params.archaic_mem.plus(task.attempt.minus(1).multiply(params.archaic_mem_ramp))+' '+params.archaic_mem_unit }
   time { params.archaic_timeout.plus(task.attempt.minus(1).multiply(params.archaic_timeout_ramp))+params.archaic_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/VCFs", mode: 'copy', pattern: '*.vcf.g{z,z.tbi}'

   input:
   tuple val(regionname), path(minmafsites), path(arcvcf), path(arctbi)

   output:
   tuple val(regionname), path("${outprefix}.vcf.gz"), path("${outprefix}.vcf.gz.tbi")

   shell:
   outprefix="${params.run_name}_archaics_${regionname}_regionspan"
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   arcregion="!{outprefix}.vcf.gz"
   echo "Extracting the MAF-filtered modern sites from the archaic VCF !{arcvcf} into ${arcregion}"
   bcftools view -T !{minmafsites} -Oz -o ${arcregion} !{arcvcf}
   tabix -f ${arcregion}
   '''
}

process vcfmerge {
   tag "${regionname}"

   cpus params.merge_cpus
   memory { params.merge_mem.plus(task.attempt.minus(1).multiply(params.merge_mem_ramp))+' '+params.merge_mem_unit }
   time { params.merge_timeout.plus(task.attempt.minus(1).multiply(params.merge_timeout_ramp))+params.merge_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/VCFs", mode: 'copy', pattern: '*.vcf.g{z,z.tbi}'

   input:
   tuple val(regionname), path(modernregionminmaf), path(moderntbi), path(arcregion), path(arctbi)

   output:
   tuple val(regionname), path("${outprefix}.vcf.gz"), path("${outprefix}.vcf.gz.tbi")

   shell:
   outprefix="${params.run_name}_merged_${regionname}_regionspan_MAFge${params.minmaf}"
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   mergedvcf="!{outprefix}.vcf.gz"
   echo "Merging the archaics (!{arcregion}) with moderns (!{modernregionminmaf}) into ${mergedvcf}"
   bcftools merge -m all -Oz -o ${mergedvcf} !{arcregion} !{modernregionminmaf}
   tabix -f ${mergedvcf}
   '''
}

process addanc {
   tag "${regionname}"

   cpus params.addanc_cpus
   memory { params.addanc_mem.plus(task.attempt.minus(1).multiply(params.addanc_mem_ramp))+' '+params.addanc_mem_unit }
   time { params.addanc_timeout.plus(task.attempt.minus(1).multiply(params.addanc_timeout_ramp))+params.addanc_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/VCFs", mode: 'copy', pattern: '*.vcf.g{z,z.tbi}'

   input:
   tuple val(regionname), path(mergedvcf), path(mergedtbi), path(ancfasta)

   output:
   tuple val(regionname), path("${outprefix}.vcf.gz"), path("${outprefix}.vcf.gz.tbi")

   shell:
   outprefix="${params.run_name}_merged_${regionname}_regionspan_MAFge${params.minmaf}_wAA"
   '''
   module load !{params.mod_samtools}
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   ancfasta="!{ancfasta}"
   renamedancfasta=$(basename ${ancfasta} .fa)_renamed.fa
   echo "Setting FASTA headers to just chromosome in ancestral allele FASTA ${ancfasta} (${renamedancfasta})"
   awk '/^>/{split(\$1, a, ":");print ">"a[3];}!/^>/{print;}' ${ancfasta} > ${renamedancfasta}
   samtools faidx ${renamedancfasta}
   aavcf="!{outprefix}.vcf.gz"
   echo "Adding ancestral allele tag for SNPs in merged VCF (!{mergedvcf}->${aavcf})"
   printf "##INFO=<ID=AA,Number=1,Type=String,Description=\\"Ancestral allele\\">\\n" > AA_header.txt
   bcftools +fill-from-fasta -Oz -o ${aavcf} !{mergedvcf} -- -i 'TYPE="SNP"' -c INFO/AA -f ${renamedancfasta} -h AA_header.txt
   tabix -f ${aavcf}
   '''
}

process hapsout {
   tag "${regionname}"

   cpus params.hapsout_cpus
   memory { params.hapsout_mem.plus(task.attempt.minus(1).multiply(params.hapsout_mem_ramp))+' '+params.hapsout_mem_unit }
   time { params.hapsout_timeout.plus(task.attempt.minus(1).multiply(params.hapsout_timeout_ramp))+params.hapsout_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/haplotypes", mode: 'copy', pattern: '*.{tsv.gz,nex}'

   input:
   tuple val(regionname), path(aavcf), path(aatbi)
   path samplemap

   output:
   tuple val(regionname), path("${outprefix}_haps_forR.tsv.gz"), emit: haptable
   tuple val(regionname), path("${outprefix}_forHapNetwork.nex"), emit: hapnexus

   shell:
   outprefix="${params.run_name}_merged_${regionname}_regionspan_MAFge${params.minmaf}_wAA"
   '''
   module load !{params.mod_bcftools}
   hapsout="!{outprefix}_haps_forR.tsv.gz"
   echo "Generating haplotypes for plotting (${hapsout})"
   bcftools query -H -f '%CHROM:%POS\\t%REF\\t%ALT\\t%INFO/AA[\\t%GT]\\n' !{aavcf} | \\
      !{projectDir}/HumanPopGenScripts/HaplotypePlotting/VCFtoRhaps.awk -v "ploidy=2" | \\
      gzip -9 > ${hapsout}
   nexout="!{outprefix}_forHapNetwork.nex"
   echo "Generating NEXUS file for haplotype network (${nexout})"
   bcftools view -v snps -c 1 -m 2 -M 2 -Ou !{aavcf} | \\
      bcftools query -H -f '%CHROM:%POS[\\t%TGT]\\n' | \\
      !{projectDir}/HumanPopGenScripts/HaplotypePlotting/query_to_haplotype_NEXUS.awk -v "idcol=!{params.samplecol}" -v "traitcol=!{params.regioncol}" - !{samplemap} > ${nexout}
   '''
}

process cadd {
   tag "${regionname}"

   cpus params.cadd_cpus
   memory { params.cadd_mem.plus(task.attempt.minus(1).multiply(params.cadd_mem_ramp))+' '+params.cadd_mem_unit }
   time { params.cadd_timeout.plus(task.attempt.minus(1).multiply(params.cadd_timeout_ramp))+params.cadd_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/annotations", mode: 'copy', pattern: '*.tsv.gz'

   input:
   tuple val(regionname), path(aavcf), path(aatbi), path(cadd)

   output:
   tuple val(regionname), path("${outprefix}.tsv.gz")

   shell:
   outprefix="${params.run_name}_merged_${regionname}_regionspan_MAFge${params.minmaf}_wAA_CADD"
   '''
   module load !{params.mod_bcftools}
   caddout="!{outprefix}.tsv.gz"
   echo "Extracting CADD scores for selected SNPs (${caddout})"
   bcftools query -f '%CHROM:%POS\\t%REF\\t%ALT\\n' !{aavcf} | \\
      !{projectDir}/HumanPopGenScripts/HaplotypePlotting/addCADD.awk - <(gzip -dc !{cadd}) | \\
      !{projectDir}/HumanPopGenScripts/HaplotypePlotting/dedupCADD.awk | \\
      gzip -9 > ${caddout}
   '''
}

process pcgs {
   tag "${regionname}"

   cpus params.pcgs_cpus
   memory { params.pcgs_mem.plus(task.attempt.minus(1).multiply(params.pcgs_mem_ramp))+' '+params.pcgs_mem_unit }
   time { params.pcgs_timeout.plus(task.attempt.minus(1).multiply(params.pcgs_timeout_ramp))+params.pcgs_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/annotations", mode: 'copy', pattern: '*.tsv.gz'

   input:
   tuple val(regionname), path(minmafsites)
   path gff

   output:
   tuple val(regionname), path("${outprefix}.tsv.gz")

   shell:
   outprefix="${params.run_name}_${regionname}_modernMAFfiltered_PCGs"
   '''
   module load !{params.mod_bedtools}
   annotout="!{outprefix}.tsv.gz"
   echo "Identifying protein-coding genes overlapping selected sites (${annotout})"
   bedtools intersect \
      -a <(!{projectDir}/HumanPopGenScripts/HaplotypePlotting/sitesToBED.awk !{minmafsites}) \
      -b <(gzip -dc !{gff} | !{projectDir}/HumanPopGenScripts/HaplotypePlotting/subsetGFF3.awk -v "feature=exon" -v "trimchr=1") \
      -wao | \\
      !{projectDir}/HumanPopGenScripts/HaplotypePlotting/overlappingFeatures.awk -v "tag=gene_name" -v "feature=exon" -v "type=protein_coding" | \\
      uniq | \\
      gzip -9 > ${annotout}
   '''
}

process vep {
   tag "${regionname}"

   cpus params.vep_cpus
   memory { params.vep_mem.plus(task.attempt.minus(1).multiply(params.vep_mem_ramp))+' '+params.vep_mem_unit }
   time { params.vep_timeout.plus(task.attempt.minus(1).multiply(params.vep_timeout_ramp))+params.vep_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/VCFs", mode: 'copy', pattern: '*.vcf'
   publishDir path: "${params.output_dir}/annotations", mode: 'copy', pattern: '*.tsv.gz'

   input:
   tuple val(regionname), path(aavcf), path(aatbi)

   output:
   tuple val(regionname), path("${outprefix}.vcf"), path("${outprefix}_worst.tsv.gz")

   shell:
   outprefix="${params.run_name}_merged_${regionname}_regionspan_MAFge${params.minmaf}_wAA_VEP"
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_vep}
   vepvcf="!{outprefix}.vcf"
   vepout="!{outprefix}_worst.tsv.gz"
   echo "Running VEP on selected sites and selecting the worst predicted consequences (${vepout})"
   cache_options="--cache --dir_cache !{params.vep_cache} --species !{params.vep_species} --assembly !{params.vep_assembly}"
   vepout_options="--vcf --regulatory --symbol --canonical --biotype --domains"
   vep ${cache_options} -i !{aavcf} -o ${vepvcf} ${vepout_options}
   bcftools +split-vep -f '%CHROM:%POS\\t%CSQ\\n' -A tab -s worst:any ${vepvcf} | \\
      gzip -9 > ${vepout}
   '''
}

process matrix {
   tag "${regionname}"

   cpus params.matrix_cpus
   memory { params.matrix_mem.plus(task.attempt.minus(1).multiply(params.matrix_mem_ramp))+' '+params.matrix_mem_unit }
   time { params.matrix_timeout.plus(task.attempt.minus(1).multiply(params.matrix_timeout_ramp))+params.matrix_timeout_unit }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/annotations", mode: 'copy', pattern: '*.tsv.gz'

   input:
   tuple val(regionname), val(chrom), val(region), val(sprimetract), val(corehap), val(pilottract), val(publishedvariants), path(sprimescore), path(corehapmatches), path(pilotscore), path(minmafsites)
   path refgenome
   path k562bed
   path jurkatbed
   path dbsnpvcf
   path dbsnptbi

   output:
   path("${params.run_name}_${regionname}_site_presence_matrix.tsv.gz")

   shell:
   '''
   module load !{params.mod_bedtools}
   module load !{params.mod_bcftools}
   printf "!{region}\\n" | awk 'BEGIN{FS="\\t";OFS=FS;}{split(\$1, a, ":");split(a[2], b, "-");print a[1], b[1]-1, b[2];}' > !{regionname}.bed
   emvars="PIB92_pilot_!{regionname}_emVars.tsv"
   echo "Extracting emVars within the region (${emvars})"
   bedtools intersect -a !{k562bed} -b !{regionname}.bed -g !{refgenome} -sorted -wa | \\
      cat - <(bedtools intersect -a !{jurkatbed} -b !{regionname}.bed -g !{refgenome} -sorted -wa) | \\
      cut -f1,3 | \\
      sort -k1,1V -k2,2n | \\
      uniq > ${emvars}
   coreseg=$(printf "!{corehap}\\n" | cut -d_ -f3 | tr -d "\n")
   core="!{regionname}_Sprime_corehap_sites.tsv"
   echo "Extracting core haplotype alleles (${core})"
   !{projectDir}/HumanPopGenScripts/HaplotypePlotting/extractTractAlleles.awk -v "chrom=!{chrom}" -v "seg=${coreseg}" !{corehapmatches} | \\
      cut -f1,2 > ${core}
   sprimeseg=$(printf "!{sprimetract}\\n" | cut -d_ -f3 | tr -d "\\n")
   sprime="!{regionname}_Sprime_sites.tsv"
   echo "Extracting Sprime tract alleles (${sprime})"
   !{projectDir}/HumanPopGenScripts/HaplotypePlotting/extractTractAlleles.awk -v "chrom=!{chrom}" -v "seg=${sprimeseg}" !{sprimescore} | \\
      cut -f1,2 > ${sprime}
   pilotseg=$(printf "!{pilottract}\\n" | cut -d_ -f3 | tr -d "\\n")
   pilot="!{regionname}_pilot_sites.tsv"
   echo "Extracting pilot Sprime tract alleles (${pilot})"
   !{projectDir}/HumanPopGenScripts/HaplotypePlotting/extractTractAlleles.awk -v "chrom=!{chrom}" -v "seg=${pilotseg}" !{pilotscore} | \\
      cut -f1,2 > ${pilot}
   sitepresencematrix="!{params.run_name}_!{regionname}_site_presence_matrix.tsv.gz"
   echo "Generating the site presence matrix (${sitepresencematrix})"
   if [[ -s "!{publishedvariants}" ]]; then
      mkdir published
      while read fn;
         do
         ln -s ${fn} published/;
      done < !{publishedvariants}
      rsidlist="!{regionname}_published_variants_rsidlist.txt"
      rsidmap="!{regionname}_published_variants_rsidmap.tsv"
      echo "Getting assembly coordinates for published SNPs based on dbSNP rsids (${rsidmap})"
      cat published/*_!{regionname}.txt | \\
         sort | \\
         uniq > ${rsidlist}
      bcftools query -i "ID=@${rsidlist}" -f '%ID\\t%CHROM\\t%POS\\n' !{dbsnpvcf} > ${rsidmap}
      corenea=$(ls published/*_!{regionname}.tsv)
      !{projectDir}/HumanPopGenScripts/HaplotypePlotting/sitepresenceMatrix.awk ${rsidmap} published/*_!{regionname}.txt ${emvars} ${corenea} ${core} ${sprime} ${pilot} !{minmafsites} | \\
         gzip -9 > ${sitepresencematrix}
   else
      !{projectDir}/HumanPopGenScripts/HaplotypePlotting/sitepresenceMatrix.awk -v "no_published=1" ${emvars} ${core} ${sprime} ${pilot} !{minmafsites} | \\
         gzip -9 > ${sitepresencematrix}
   fi
   '''
}

workflow {
   //Load the samples to keep file:
   samples_to_keep = file(params.samples_to_keep, checkIfExists: true)

   //Load the sample to region map:
   samplemap = file(params.samplemap, checkIfExists: true)

   //Load the annotation GFF3:
   ref_annotation = file(params.ref_annotation, checkIfExists: true)

   //Load the dbSNP VCF and index:
   dbsnp_vcf = file(params.dbsnp_vcf, checkIfExists: true)
   dbsnp_tbi = file(params.dbsnp_vcf+'.tbi', checkIfExists: true)

   //Load the ref .genome file for faster BEDtools calls:
   refgenome = file(params.ref.replaceFirst("[.]fn?a(sta)?([.]gz)?", ".genome"), checkIfExists: true)

   //Load the MPRA emVar BEDs:
   k562_bed = file(params.k562_bed, checkIfExists: true)
   jurkat_bed = file(params.jurkat_bed, checkIfExists: true)

   //Load the target regions:
   Channel
      .fromPath(params.target_regions, checkIfExists: true)
      .ifEmpty { error "Unable to find target regions file: ${params.target_regions}" }
      .splitCsv(sep: "\t", header: false)
      .tap { target_regions }
      .map { [ it[0], it[1], it[2] ] }
      .tap { target_regions_min }
      .subscribe { println "Added region ${it[1]} to target_regions channel" }

   //Load the per-chromosome modern sample VCFs:
   Channel
      .fromPath(params.modern_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find modern sample VCFs matching glob: ${params.modern_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.modern_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from modern sample VCF: ${params.modern_regex}" }
      .tap { modern_vcfs_only }
      .subscribe { println "Added ${it[1]} to modern_vcfs_only channel" }
   Channel
      .fromPath(params.modern_glob+'.tbi', checkIfExists: true)
      .ifEmpty { error "Unable to find modern sample VCF indices matching glob: ${params.modern_glob}.tbi" }
      .map { a -> [ (a.getSimpleName() =~ params.modern_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from modern sample VCF index: ${params.modern_regex}" }
      .tap { modern_tbis_only }
      .subscribe { println "Added ${it[1]} to modern_tbis_only channel" }
   modern_vcfs = modern_vcfs_only
      .join(modern_tbis_only, by: 0, failOnDuplicate: true, failOnMismatch: true)

   //Load the per-chromosome archaic hominin VCFs:
   Channel
      .fromPath(params.archaic_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find archaic sample VCFs matching glob: ${params.archaic_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.archaic_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from archaic sample VCF: ${params.archaic_regex}" }
      .tap { archaic_vcfs_only }
      .subscribe { println "Added ${it[1]} to archaic_vcfs_only channel" }
   Channel
      .fromPath(params.archaic_glob+'.tbi', checkIfExists: true)
      .ifEmpty { error "Unable to find archaic sample VCF indices matching glob: ${params.archaic_glob}.tbi" }
      .map { a -> [ (a.getSimpleName() =~ params.archaic_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from archaic sample VCF index: ${params.archaic_regex}" }
      .tap { archaic_tbis_only }
      .subscribe { println "Added ${it[1]} to archaic_tbis_only channel" }
   archaic_vcfs = archaic_vcfs_only
      .join(archaic_tbis_only, by: 0, failOnDuplicate: true, failOnMismatch: true)

   //Load the ancestral allele FASTAs:
   Channel
      .fromPath(params.aa_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find ancestral allele FASTAs matching glob: ${params.aa_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.aa_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from ancestral allele FASTA: ${params.aa_regex}" }
      .tap { anc_fastas }
      .subscribe { println "Added ${it[1]} to anc_fastas channel" }

   //Load the CADD TSVs from WGSA:
   Channel
      .fromPath(params.cadd_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find WGSA CADD TSVs matching glob: ${params.cadd_glob}" }
      .map { a -> [ (a.getName() =~ params.cadd_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from WGSA CADD TSV: ${params.cadd_regex}" }
      .tap { cadd_tsvs }
      .subscribe { println "Added ${it[1]} to cadd_tsvs channel" }

   //Load the final Sprime score files:
   //Regex capture groups: 1=Population, 2=Chromosome
   Channel
      .fromPath(params.final_sprime_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find final Sprime score files matching glob: ${params.final_sprime_glob}" }
      .map { a -> [ (a.getBaseName() =~ params.final_sprime_regex)[0][2], (a.getBaseName() =~ params.final_sprime_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID and population from final Sprime score file: ${params.final_sprime_regex}" }
      .tap { final_sprime_scores }
      .subscribe { println "Added ${it[2]} to final_sprime_scores channel" }

   //Load the core haplotype matches files:
   //Regex capture groups: 1=Population
   Channel
      .fromPath(params.corehap_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find corehap matches files matching glob: ${params.corehap_glob}" }
      .map { a -> [ (a.getBaseName() =~ params.corehap_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract population from corehap matches file: ${params.corehap_regex}" }
      .tap { corehap_matches }
      .subscribe { println "Added ${it[1]} to corehap_matches channel" }

   //Load the pilot Sprime score files:
   //Regex capture groups: 1=Chromosome, 2=Population
   Channel
      .fromPath(params.pilot_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find pilot Sprime score files matching glob: ${params.pilot_glob}" }
      .map { a -> [ (a.getBaseName() =~ params.pilot_regex)[0][1], (a.getBaseName() =~ params.pilot_regex)[0][2], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID and population from pilot Sprime score file: ${params.pilot_regex}" }
      .tap { pilot_sprime_scores }
      .subscribe { println "Added ${it[2]} to pilot_sprime_scores channel" }

   //Extract variants from modern samples in the region and apply MAF filter:
   modern_inputs = target_regions_min
      .combine(modern_vcfs, by: 0)
      .map({ [ it[1], it[2], it[3], it[4] ] })
   modern(modern_inputs, samples_to_keep)

   //Extract the same variants from the archaics:
   archaic_inputs = modern.out.mafsites
      .join(target_regions_min
         .combine(archaic_vcfs, by: 0)
         .map({ [ it[1], it[3], it[4] ] }), by: 0, failOnDuplicate: true, failOnMismatch: true)
   archaic(archaic_inputs)

   //Merge the modern and archaic genotypes:
   merge_inputs = modern.out.mafvcf
      .join(archaic.out, by: 0, failOnDuplicate: true, failOnMismatch: true)
   vcfmerge(merge_inputs)

   //Add ancestral allele calls as INFO/AA to the merged VCF:
   anc_inputs = vcfmerge.out
      .join(target_regions_min
         .combine(anc_fastas, by: 0)
         .map({ [ it[1], it[3] ] }), by: 0, failOnDuplicate: true, failOnMismatch: true)
   addanc(anc_inputs)

   //Generate the haplotype outputs for plotting and haplotype networks:
   hapsout(addanc.out, samplemap)

   //Extract CADD score annotations for the variants:
   cadd_inputs = addanc.out
      .join(target_regions_min
         .combine(cadd_tsvs, by: 0)
         .map({ [ it[1], it[3] ] }), by: 0, failOnDuplicate: true, failOnMismatch: true)
   cadd(cadd_inputs)

   //Annotate the variants that overlap with exons of protein-coding genes:
   pcgs(modern.out.mafsites, ref_annotation)

   //Perform VEP annotation of the variants:
   vep(addanc.out)

   //Construct the site-presence matrix:
   final_score = target_regions
      .map({ [ it[0], (it[3] =~ /^([0-9A-Za-z-,]+)_(\p{Alnum}+)_([0-9]+)$/)[0][1], it[1] ] })
      .combine(final_sprime_scores, by: [0,1])
      .map({ [ it[2], it[3] ] })
   corehap_score = target_regions
      .map({ [ (it[4] =~ /^([0-9A-Za-z-,]+)_(\p{Alnum}+)_([0-9A-Za-z.]+)$/)[0][1], it[1] ] })
      .combine(corehap_matches, by: 0)
      .map({ [ it[1], it[2] ] })
   pilot_score = target_regions
      .map({ [ it[0], (it[5] =~ /^([0-9A-Za-z-]+)_(\p{Alnum}+)_([0-9]+)$/)[0][1], it[1] ] })
      .combine(pilot_sprime_scores, by: [0,1])
      .map({ [ it[2], it[3] ] })
   matrix_inputs = target_regions
      .map({ [ it[1], it[0], it[2], it[3], it[4], it[5], it[6] ] })
      .join(final_score, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .join(corehap_score, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .join(pilot_score, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .join(modern.out.mafsites, by: 0, failOnDuplicate: true, failOnMismatch: true)
   matrix(matrix_inputs, refgenome, k562_bed, jurkat_bed, dbsnp_vcf, dbsnp_tbi)
}

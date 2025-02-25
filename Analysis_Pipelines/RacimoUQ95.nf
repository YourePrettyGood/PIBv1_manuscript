#!/usr/bin/env nextflow
/* Pipeline to evaluate the Racimo et al. 2017 U and Q95 statistics for *
 *  adaptive introgression.                                             *
 * Steps:                                                               *
 *  Reformat ancestral allele FASTA headers to match ref chroms ->      *
 *  Extract sites where archaic genotypes are not missing ->            *
 *  Identify sites where archaics match y == 1.0 || z == 1.0 ->         *
 *  Combine AMB, DEN, and NEA sites together for filtering moderns ->   *
 *  Extract the above SNPs from moderns and re-calculate per-pop AFs -> *
 *  Construct the table of w, x, y, and z for each SNP ->               *
 *  Evaluate U and Q95 in windows based on the above table              */

nextflow.enable.dsl=2

//Default paths, globs, and regexes:
//Modern sample VCF glob:
params.modern_vcf_glob = "${projectDir}/modern_VCFs/*.vcf.gz"
//Regex to extract chromosome number from modern sample VCF:
params.modern_vcf_regex = ~/_chr(\p{Alnum}+)/
//Archaic sample VCF glob:
params.arc_vcf_glob = "${projectDir}/archaic_VCFs/*.vcf.gz"
//Regex to extract chromosome number from archaic sample VCF:
params.arc_vcf_regex = ~/_chr(\p{Alnum}+)/

//U and Q95 statistic window size:
params.window_size = "40000"
//Use population/analsis group labels instead of regions?:
params.usepop = "0"
//Populations to exclude from allele frequency calculations:
//e.g. ASW for 1000 Genomes, perhaps ASW,ACB as well
params.pops_to_exclude = ""

//Defaults for cpus, memory, time, and queue for each process:
//AA FASTA reformat (aachr):
params.aachr_cpus = 1
params.aachr_mem = 1
params.aachr_mem_increment = 4
params.aachr_timeout = '1h'
params.aachr_retry_timeout = '4h'
params.aachr_queue = 'day'
params.aachr_retry_queue = 'day'
//Extract archaic present sites (arcpresent):
params.arcpresent_cpus = 1
params.arcpresent_mem = 1
params.arcpresent_mem_increment = 4
params.arcpresent_timeout = '1h'
params.arcpresent_retry_timeout = '4h'
params.arcpresent_queue = 'day'
params.arcpresent_retry_queue = 'day'
//Extract archaic derived sites (yzone):
params.yzone_cpus = 1
params.yzone_mem = 1
params.yzone_mem_increment = 4
params.yzone_timeout = '2h'
params.yzone_retry_timeout = '6h'
params.yzone_queue = 'day'
params.yzone_retry_queue = 'day'
//Concatenate sites across archaic origins (arcconcat):
params.arcconcat_cpus = 1
params.arcconcat_mem = 1
params.arcconcat_mem_increment = 4
params.arcconcat_timeout = '1h'
params.arcconcat_retry_timeout = '4h'
params.arcconcat_queue = 'day'
params.arcconcat_retry_queue = 'day'
//Extract modern genotypes at archaic derived sites (extractyzone):
params.extractyzone_cpus = 1
params.extractyzone_mem = 1
params.extractyzone_mem_increment = 4
params.extractyzone_timeout = '1h'
params.extractyzone_retry_timeout = '4h'
params.extractyzone_queue = 'day'
params.extractyzone_retry_queue = 'day'
//Construct table of modern and archaic allele frequencies (joinwxyz):
params.joinwxyz_cpus = 1
params.joinwxyz_mem = 4
params.joinwxyz_mem_increment = 4
params.joinwxyz_timeout = '1h'
params.joinwxyz_retry_timeout = '4h'
params.joinwxyz_queue = 'day'
params.joinwxyz_retry_queue = 'day'
//Evaluate U and Q95 statistics (evaluq):
params.evaluq_cpus = 1
params.evaluq_mem = 1
params.evaluq_mem_increment = 4
params.evaluq_timeout = '1h'
params.evaluq_retry_timeout = '4h'
params.evaluq_queue = 'day'
params.evaluq_retry_queue = 'day'
//Concatenate the U and Q95 statistic outputs (cat_outs):
params.cat_outs_cpus = 1
params.cat_outs_mem = 2
params.cat_outs_mem_increment = 4
params.cat_outs_timeout = '1h'
params.cat_outs_retry_timeout = '4h'
params.cat_outs_queue = 'day'
params.cat_outs_retry_queue = 'day'

process aachr {
   tag "${chrom}"

   cpus params.aachr_cpus
   memory { params.aachr_mem.plus(task.attempt.minus(1).multiply(params.aachr_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.aachr_retry_timeout : params.aachr_timeout }
   queue { task.attempt >= 2 ? params.aachr_retry_queue : params.aachr_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(ancfasta)

   output:
   path("AArename_${chrom}.stderr"), emit: logs
   tuple val(chrom), path("anc_reheadered_${chrom}.fa"), emit: fa

   shell:
   '''
   !{projectDir}/HumanPopGenScripts/AdaptiveIntrogression/fixAAFASTAname.awk !{ancfasta} > anc_reheadered_!{chrom}.fa 2> AArename_!{chrom}.stderr
   '''
}

process arcpresent {
   tag "${chrom}"

   cpus params.arcpresent_cpus
   memory { params.arcpresent_mem.plus(task.attempt.minus(1).multiply(params.arcpresent_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.arcpresent_retry_timeout : params.arcpresent_timeout }
   queue { task.attempt >= 2 ? params.arcpresent_retry_queue : params.arcpresent_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(arcvcf), path(arctbi), path(ancfasta)
   path(arc_pop_map)
   path(aa_header)

   output:
   tuple path("bcftools_filltags_arc_${params.run_name}_${chrom}.stderr"), path("bcftools_annotate_trimarctags_${params.run_name}_${chrom}.stderr"), path("bcftools_fillfromfasta_archaics_${params.run_name}_${chrom}.stderr"), emit: logs
   tuple val(chrom), path("${params.run_name}_ArchaicsPresent_AA_${chrom}.vcf.gz"), path("${params.run_name}_ArchaicsPresent_AA_${chrom}.vcf.gz.tbi"), emit: vcf

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   bcftools +fill-tags -Ou !{arcvcf} -- -t AC,AF,AN -S !{arc_pop_map} 2> bcftools_filltags_arc_!{params.run_name}_!{chrom}.stderr | \\
      bcftools annotate -Ou -i "AN_DEN>0&&AN_NEA>0" -x ^FORMAT/GT,FORMAT/GQ 2> bcftools_annotate_trimarctags_!{params.run_name}_!{chrom}.stderr | \\
      bcftools +fill-from-fasta -Oz -o !{params.run_name}_ArchaicsPresent_AA_!{chrom}.vcf.gz - -- -c AA -f !{ancfasta} -h !{aa_header} 2> bcftools_fillfromfasta_archaics_!{params.run_name}_!{chrom}.stderr
   tabix -f !{params.run_name}_ArchaicsPresent_AA_!{chrom}.vcf.gz
   '''
}

process yzone {
   tag "${chrom} ${origin}"

   cpus params.yzone_cpus
   memory { params.yzone_mem.plus(task.attempt.minus(1).multiply(params.yzone_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.yzone_retry_timeout : params.yzone_timeout }
   queue { task.attempt >= 2 ? params.yzone_retry_queue : params.yzone_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(origin), val(chrom), path(arcpresentvcf), path(arcpresenttbi)

   output:
   path("keepyz1_${origin}_${chrom}.stderr"), emit: logs
   tuple val(chrom), path("${params.run_name}_${origin}_${chrom}.vcf.gz"), path("${params.run_name}_${origin}_${chrom}.vcf.gz.tbi"), emit: vcf

   shell:
   '''
   module load !{params.mod_htslib}
   bgzip -dc !{arcpresentvcf} | \\
      !{projectDir}/HumanPopGenScripts/AdaptiveIntrogression/RacimoUQ95_keepyz1.awk -v "origin=!{origin}" 2> keepyz1_!{origin}_!{chrom}.stderr | \\
      bgzip -c > !{params.run_name}_!{origin}_!{chrom}.vcf.gz
   tabix -f !{params.run_name}_!{origin}_!{chrom}.vcf.gz
   '''
}

process arcconcat {
   tag "${chrom}"

   cpus params.arcconcat_cpus
   memory { params.arcconcat_mem.plus(task.attempt.minus(1).multiply(params.arcconcat_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.arcconcat_retry_timeout : params.arcconcat_timeout }
   queue { task.attempt >= 2 ? params.arcconcat_retry_queue : params.arcconcat_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(arcvcfs), path(arctbis)

   output:
   path("bcftools_concat_${params.run_name}_${chrom}.stderr"), emit: logs
   tuple val(chrom), path("${params.run_name}_archaics_AMBDENNEA_${chrom}.vcf.gz"), path("${params.run_name}_archaics_AMBDENNEA_${chrom}.vcf.gz.tbi"), emit: vcf

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   bcftools concat -a -Oz -o !{params.run_name}_archaics_AMBDENNEA_!{chrom}.vcf.gz !{arcvcfs} 2> bcftools_concat_!{params.run_name}_!{chrom}.stderr
   tabix -f !{params.run_name}_archaics_AMBDENNEA_!{chrom}.vcf.gz
   '''
}

process extractyzone {
   tag "${chrom}"

   cpus params.extractyzone_cpus
   memory { params.extractyzone_mem.plus(task.attempt.minus(1).multiply(params.extractyzone_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.extractyzone_retry_timeout : params.extractyzone_timeout }
   queue { task.attempt >= 2 ? params.extractyzone_retry_queue : params.extractyzone_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(arcvcf), path(arctbi), path(modernvcf), path(moderntbi)
   path(modern_pop_map)

   output:
   tuple path("bcftools_view_filtermodern_${params.run_name}_${chrom}.stderr"), path("bcftools_annotate_trimmoderntags_${params.run_name}_${chrom}.stderr"), path("bcftools_filltags_modern_${params.run_name}_${chrom}.stderr"), emit: logs
   tuple val(chrom), path("${params.run_name}_yz1_ACAFAN_${chrom}.vcf.gz"), path("${params.run_name}_yz1_ACAFAN_${chrom}.vcf.gz.tbi"), emit: vcf

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   bcftools view -Ou -v snps -T !{arcvcf} !{modernvcf} 2> bcftools_view_filtermodern_!{params.run_name}_!{chrom}.stderr | \\
      bcftools annotate -Ou -x INFO,FORMAT - 2> bcftools_annotate_trimmoderntags_!{params.run_name}_!{chrom}.stderr | \\
      bcftools +fill-tags -Oz -o !{params.run_name}_yz1_ACAFAN_!{chrom}.vcf.gz - -- -t AC,AF,AN -S !{modern_pop_map} 2> bcftools_filltags_modern_!{params.run_name}_!{chrom}.stderr
   tabix -f !{params.run_name}_yz1_ACAFAN_!{chrom}.vcf.gz
   '''
}

process joinwxyz {
   tag "${chrom}"

   cpus params.joinwxyz_cpus
   memory { params.joinwxyz_mem.plus(task.attempt.minus(1).multiply(params.joinwxyz_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.joinwxyz_retry_timeout : params.joinwxyz_timeout }
   queue { task.attempt >= 2 ? params.joinwxyz_retry_queue : params.joinwxyz_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(arcvcf), path(arctbi), path(modernvcf), path(moderntbi)

   output:
   path("joinwxyz_${params.run_name}_${chrom}.stderr"), emit: logs
   tuple val(chrom), path("${params.run_name}_UQ95_siteAFs_${chrom}.tsv.gz"), emit: wxyz

   shell:
   '''
   module load !{params.mod_htslib}
   !{projectDir}/HumanPopGenScripts/AdaptiveIntrogression/RacimoUQ95_joinwxyz.awk -v "usepop=!{params.usepop}" -v "exclude=!{params.pops_to_exclude}" <(bgzip -dc !{arcvcf}) <(bgzip -dc !{modernvcf}) 2> joinwxyz_!{params.run_name}_!{chrom}.stderr | \\
      gzip -9 > !{params.run_name}_UQ95_siteAFs_!{chrom}.tsv.gz
   '''
}

process evaluq {
   tag "${chrom} ${A} ${B} ${origin}"

   cpus params.evaluq_cpus
   memory { params.evaluq_mem.plus(task.attempt.minus(1).multiply(params.evaluq_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.evaluq_retry_timeout : params.evaluq_timeout }
   queue { task.attempt >= 2 ? params.evaluq_retry_queue : params.evaluq_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(chrom), path(wxyz), val(A), val(B), val(C), val(D), val(w), val(x), val(origin)

   output:
   path("UQ95_calculate_${params.run_name}_w${params.window_size}_${A}_${B}_${origin}_w${w}_x${x}_${chrom}.stderr"), emit: logs
   path("${params.run_name}_UQ95_w${params.window_size}_${A}_${B}_${origin}_w${w}_x${x}_${chrom}.tsv"), emit: tsv

   shell:
   '''
   if [[ "!{origin}" == "AMB" ]]; then
      y="1.0";
      z="1.0";
   elif [[ "!{origin}" == "DEN" ]]; then
      y="0.0";
      z="1.0";
   elif [[ "!{origin}" == "NEA" ]]; then
      y="1.0";
      z="0.0";
   else
      echo "Error: Archaic origin label !{origin} not in AMB,DEN,NEA, cannot proceed." 1>&2;
      exit 2;
   fi
   gzip -dc !{wxyz} | \\
      !{projectDir}/HumanPopGenScripts/AdaptiveIntrogression/RacimoUQ95_calculate.awk -v "A=!{A}" -v "B=!{B}" -v "C=!{C}" -v "D=!{D}" -v "w=!{w}" -v "x=!{x}" -v "y=${y}" -v "z=${z}" -v "arclabel=!{origin}" -v "window_size=!{params.window_size}" -v "header=1" 2> UQ95_calculate_!{params.run_name}_w!{params.window_size}_!{A}_!{B}_!{origin}_w!{w}_x!{x}_!{chrom}.stderr > !{params.run_name}_UQ95_w!{params.window_size}_!{A}_!{B}_!{origin}_w!{w}_x!{x}_!{chrom}.tsv
   '''
}

process cat_outs {
   cpus params.cat_outs_cpus
   memory { params.cat_outs_mem.plus(task.attempt.minus(1).multiply(params.cat_outs_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.cat_outs_retry_timeout : params.cat_outs_timeout }
   queue { task.attempt >= 2 ? params.cat_outs_retry_queue : params.cat_outs_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/UQ95", mode: 'copy', pattern: '*.tsv'

   input:
   path(tsvs)

   output:
   path("${params.run_name}_UQ95_w${params.window_size}.tsv"), emit: tsv

   shell:
   '''
   #Symlink the many input files to avoid SLURM SBATCH script size limits:
   #Also concatenate the files without extra/redundant headers:
   header="1";
   while IFS=$'\\t' read -a a;
      do
      ln -s ${a[2]} ${a[1]};
      if [[ "${header}" == "1" ]]; then
         cat ${a[1]};
         header="0";
      else
         tail -n+2 ${a[1]};
      fi;
   done < <(sort -k2,2V < !{tsvs}) > !{params.run_name}_UQ95_w!{params.window_size}.tsv
   '''
}

workflow {
   //Detection mechanism for "chr" prefixes of autosomes:
   has_chr_prefix = params.autosomes.count('chr') > 0
   autosome_list = params.autosomes
   autosome_num_list = has_chr_prefix ? params.autosomes.replaceAll('chr', '') : autosome_list
   num_autosomes = autosome_num_list.tokenize(',').size()

   //Set up a fixed channel of archaic origins for use with yzone:
   arc_origin_list = [ "AMB", "DEN", "NEA" ]
   arc_origins = Channel
      .fromList(arc_origin_list)
      .tap { arc_origins_filter }

   //Set up the channel of U and Q95 scans to run:
   //i.e. parameter sets (A, B, C, D, w, x, origin)
   Channel
      .fromPath(params.UQ95_scan_params, checkIfExists: true)
      .ifEmpty { error "Unable to find U/Q95 scan parameters file: ${params.UQ95_scan_params}" }
      .splitCsv(sep: "\t")
      .tap { scan_params }
      .subscribe { println "Added ${it} to scan_params channel" }

   //Set up the channel of modern VCFs and indices:
   Channel
      .fromPath(params.modern_vcf_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find VCFs of modern samples matching glob: ${params.modern_vcf_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.modern_vcf_regex)[0][1], a ] }
      .ifEmpty { error "Regex failed to extract chromosome ID from modern sample VCF: ${params.modern_vcf_glob} =~ ${params.modern_vcf_regex}" }
      .tap { modern_vcfs_only }
      .subscribe { println "Added ${it[1]} to modern_vcfs_only channel" }
   Channel
      .fromPath(params.modern_vcf_glob+'.tbi', checkIfExists: true)
      .ifEmpty { error "Unable to find VCF indices of modern samples matching glob: ${params.modern_vcf_glob}.tbi" }
      .map { a -> [ (a.getSimpleName() =~ params.modern_vcf_regex)[0][1], a ] }
      .ifEmpty { error "Regex failed to extract chromosome ID from modern sample VCF index: ${params.modern_vcf_glob}.tbi =~ ${params.modern_vcf_regex}" }
      .tap { modern_tbis_only }
      .subscribe { println "Added ${it[1]} to modern_tbis_only channel" }
   modern_vcfs = modern_vcfs_only
      .join(modern_tbis_only, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })

   //Set up the channel of archaic VCFs and indices:
   Channel
      .fromPath(params.arc_vcf_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find VCFs of archaic samples matching glob: ${params.arc_vcf_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.arc_vcf_regex)[0][1], a ] }
      .ifEmpty { error "Regex failed to extract chromosome ID from archaic sample VCF: ${params.arc_vcf_glob} =~ ${params.arc_vcf_regex}" }
      .tap { arc_vcfs_only }
      .subscribe { println "Added ${it[1]} to arc_vcfs_only channel" }
   Channel
      .fromPath(params.arc_vcf_glob+'.tbi', checkIfExists: true)
      .ifEmpty { error "Unable to find VCF indices of archaic samples matching glob: ${params.arc_vcf_glob}.tbi" }
      .map { a -> [ (a.getSimpleName() =~ params.arc_vcf_regex)[0][1], a ] }
      .ifEmpty { error "Regex failed to extract chromosome ID from archaic sample VCF index: ${params.arc_vcf_glob}.tbi =~ ${params.arc_vcf_regex}" }
      .tap { arc_tbis_only }
      .subscribe { println "Added ${it[1]} to arc_tbis_only channel" }
   arc_vcfs = arc_vcfs_only
      .join(arc_tbis_only, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })

   //Set up the channel of ancestral allele FASTAs:
   Channel
      .fromPath(params.anc_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find ancestral allele FASTAs matching glob: ${params.anc_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.anc_regex)[0][1], a ] }
      .ifEmpty { error "Regex failed to extract chromosome ID from ancestral allele FASTA: ${params.anc_glob} =~ ${params.anc_regex}" }
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
      .tap { ancestral }
      .subscribe { println "Added chr${it[0]} ancestral allele FASTA (${it[1]}) to ancestral channel" }

   //Set up the value channels of the pop maps (archaic and modern):
   modern_pop_map = file(params.modern_pop_map, checkIfExists: true)
   arc_pop_map = file(params.arc_pop_map, checkIfExists: true)

   //Set up the value channel for the AA tag header:
   aa_header = file(params.aa_header, checkIfExists: true)

   //Reformat ancestral allele FASTA headers:
   aachr(ancestral)

   //Extract sites where the archaics have genotypes:
   arcpresent_inputs = arc_vcfs
      .join(aachr.out.fa, by: 0, failOnDuplicate: true, failOnMismatch: true)
   arcpresent(arcpresent_inputs, arc_pop_map, aa_header)

   //Identify sites where y == 1.0 || z == 1.0:
   yzone_inputs = arc_origins
      .combine(arcpresent.out.vcf)
   yzone(yzone_inputs)

   //Combine AMB, DEN, and NEA yz1 sites:
   arcconcat_inputs = yzone.out.vcf
      .groupTuple(by: 0, sort: true, size: arc_origin_list.size())
   arcconcat(arcconcat_inputs)

   //Extract yz1 sites from modern VCFs:
   extractyzone_inputs = arcconcat.out.vcf
      .tap { arc_yzone_vcfs }
      .join(modern_vcfs, by: 0, failOnDuplicate: true, failOnMismatch: true)
   extractyzone(extractyzone_inputs, modern_pop_map)

   //Construct the table of w, x, y, and z per SNP:
   joinwxyz_inputs = arc_yzone_vcfs
      .join(extractyzone.out.vcf, by: 0, failOnDuplicate: true, failOnMismatch: true)
   joinwxyz(joinwxyz_inputs)

   //Calculate U and Q95 statistics in windows:
   evaluq_inputs = joinwxyz.out.wxyz
      .combine(scan_params)
   evaluq(evaluq_inputs)

   //Concatenate all the results TSVs together, skipping extra headers:
   cat_outs_inputs = evaluq.out.tsv
      .collectFile() { [ "UQ95_results_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   cat_outs(cat_outs_inputs)
}

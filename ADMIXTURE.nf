#!/usr/bin/env nextflow
/* Pipeline for parallelizing ADMIXTURE across K and replicates          *
 * Steps:                                                                *
 *  Run ADMIXTURE per K and per replicate ->                             *
 *  Zip outputs per K, then zip the zips together ->                     *
 *  Run CLUMPAK in a very particular way                                 */

nextflow.enable.dsl = 2

//Default parameters:
params.minK = 2
params.maxK = 12
params.replicates = 20
params.cv_folds = 20
params.mclthreshold = 0.9

//Default paths, globs, and regexes:
//Input PLINK files:
params.plink_prefix = ""
//Input CLUMPAK install zip:
params.clumpak_zip = ""
//Path within the zip to the CLUMPAK Perl files:
params.clumpak_path = "26_03_2015_CLUMPAK/CLUMPAK"

//Defaults for cpus, memory, and time for each process:
//PLINK subsetting step:
params.subset_cpus = 4
params.subset_mem = 3
params.subset_timeout = '1h'
params.subset_queue = 'day'
//ADMIXPIPE ADMIXTURE step:
params.admixture_cpus = 1
params.admixture_mem = 16
params.admixture_timeout = '24h'
params.admixture_long_timeout = '168h'
params.admixture_queue = 'day'
params.admixture_long_queue = 'week'
//Zipping of Q files for CLUMPAK:
params.zip_q_cpus = 1
params.zip_q_mem = 1
params.zip_q_timeout = '1h'
params.zip_q_queue = 'day'
//CLUMPAK:
params.clumpak_cpus = 1
params.clumpak_mem = 2
params.clumpak_timeout = '3h'
params.clumpak_queue = 'day'

process subset {
   cpus params.subset_cpus
   memory { params.subset_mem.plus(1).plus(task.attempt.minus(1).multiply(16))+' GB' }
   time { task.attempt >= 2 ? '24h' : params.subset_timeout }
   queue params.subset_queue
   errorStrategy { task.exitStatus in ([1]+(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   path plink_bed
   path plink_bim
   path plink_fam
   path pop_map
   path excludesamples

   output:
   tuple path("${params.plink_prefix}_plink_subset.stderr"), path("${params.plink_prefix}_plink_subset.stdout"), emit: logs
   tuple path("${params.plink_prefix}_noexcluded.bed"), path("${params.plink_prefix}_noexcluded.bim"), path("${params.plink_prefix}_noexcluded.fam"), emit: plink
   path("${params.plink_prefix}_noexcluded.fam"), emit: fam
   path("admixture_pop_map.txt"), emit: pop_map

   shell:
   plink_mem = task.memory.minus(1.GB).toMega()
   plink_prefix = plink_bed.getBaseName()
   '''
   module load !{params.mod_plink}
   #Filter the population map for ADMIXTURE:
   !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excludesamples} !{pop_map} > admixture_pop_map.txt
   #Filter the samples from the PLINK files:
   !{projectDir}/HumanPopGenScripts/excludeSamples.awk -v "samplecol=2" -v "negate=1" !{excludesamples} !{plink_fam} | \
      cut -f1,2 > samples_to_exclude_forPLINK.tsv
   plink --threads !{task.cpus} --memory !{plink_mem} --bfile !{plink_prefix} --remove samples_to_exclude_forPLINK.tsv --make-bed --out !{params.plink_prefix}_noexcluded 2> !{params.plink_prefix}_plink_subset.stderr > !{params.plink_prefix}_plink_subset.stdout
   '''
}

process admixture {
   tag "K${K} r${replicate}"

   cpus params.admixture_cpus
   memory { params.admixture_mem.plus(task.attempt.minus(1).multiply(32))+' GB' }
   time { task.attempt >= 2 ? params.admixture_long_timeout : params.admixture_timeout }
   queue { task.attempt >= 2 ? params.admixture_long_queue : params.admixture_queue }
   errorStrategy { task.exitStatus in ([1]+(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/ADMIXTURE", mode: 'copy', pattern: '*.{P,Q}'

   input:
   tuple path(plink_bed), path(plink_bim), path(plink_fam), val(K), val(replicate), val(prng_seed)

   output:
   tuple path("${params.plink_prefix}_admixture_K${K}_r${replicate}.stderr"), path("${params.plink_prefix}_admixture_K${K}_r${replicate}.stdout"), emit: logs
   path("${params.plink_prefix}.${K}_${replicate}.P"), emit: P
   path("${params.plink_prefix}.${K}_${replicate}.Q"), emit: Q
   path("${params.plink_prefix}_admixture_K${K}_r${replicate}.stdout"), emit: CV

   shell:
   plink_prefix = plink_bed.getBaseName()
   '''
   module load !{params.mod_admixture}
   admixture -j!{task.cpus} --seed=!{prng_seed} --cv=!{params.cv_folds} !{plink_prefix}.bed !{K} 2> !{params.plink_prefix}_admixture_K!{K}_r!{replicate}.stderr > !{params.plink_prefix}_admixture_K!{K}_r!{replicate}.stdout
   mv !{plink_prefix}.!{K}.P !{params.plink_prefix}.!{K}_!{replicate}.P
   mv !{plink_prefix}.!{K}.Q !{params.plink_prefix}.!{K}_!{replicate}.Q
   '''
}

process zip_Q {
   cpus params.zip_q_cpus
   memory { params.zip_q_mem.plus(task.attempt.minus(1).multiply(1))+' GB' }
   time { task.attempt >= 2 ? '12h' : params.zip_q_timeout }
   queue params.zip_q_queue
   errorStrategy { task.exitStatus in ([1]..(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/ADMIXTURE", mode: 'copy', pattern: '*.zip'

   input:
   path("admixture_q_paths.tsv")
   path(plink_fam)
   path(pop_map)

   output:
   path("${params.plink_prefix}.zip"), emit: zip

   shell:
   '''
   while read -a a;
      do
      !{projectDir}/HumanPopGenScripts/ADMIXTURE/rearrangeQ.awk !{plink_fam} !{pop_map} ${a[2]} > ${a[1]};
   done < admixture_q_paths.tsv
   for K in {!{params.minK}..!{params.maxK}};
      do
      zip !{params.plink_prefix}_K${K}.zip $(find . -name "*.Q" -print | awk -v "K=${K}" '{n=split($1, a, ".");m=split(a[n-1], b, "_");if (b[1] == K) {printf " %s", $1;};}')
   done
   zip !{params.plink_prefix}.zip !{params.plink_prefix}_K*.zip
   '''
}

process clumpak {
   cpus params.clumpak_cpus
   memory { params.clumpak_mem.plus(task.attempt.minus(1).multiply(1))+' GB' }
   time { task.attempt == 2 ? '12h' : params.clumpak_timeout }
   queue params.clumpak_queue
   errorStrategy { task.exitStatus in ([1]..(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/ADMIXTURE", mode: 'copy', pattern: '*.tar.gz'

   input:
   path(CLUMPAK)
   path(pop_map)
   path(Qzip)
   path(drawparams)

   output:
   tuple path("${params.clumpak_path}/CLUMPAK_${params.plink_prefix}.stderr"), path("${params.clumpak_path}/CLUMPAK_${params.plink_prefix}.stdout"), emit: logs
   path("CLUMPAK_${params.plink_prefix}_output.tar.gz"), emit: targz

   shell:
   customdraw = drawparams.exists() ? "--drawparams ${drawparams}" : ""
   '''
   #Load the CLUMPAK module that loads the appropriate Perl and sets PERL5LIB appropriately:
   module load !{params.mod_clumpak}
   #Stage the CLUMPAK installation:
   unzip !{CLUMPAK}
   chmod u+x !{params.clumpak_path}/CLUMPP/CLUMPP
   chmod u+x !{params.clumpak_path}/distruct/distruct1.1
   chmod u+x !{params.clumpak_path}/mcl/bin/*
   #Set up the input directory and files:
   mkdir input
   cut -f2 !{pop_map} > input/clumpak_pop_map.txt
   ln -s ../!{Qzip} input/
   if [[ -n "!{customdraw}" ]]; then
      cp !{drawparams} !{params.clumpak_path}/
   fi
   #Set up the output directory:
   mkdir output
   #Capture the absolute path of the base directory so we can indicate absolute paths in the call:
   base_dir=$(pwd)
   #Set the CWD to the staged CLUMPAK install:
   pushd !{params.clumpak_path}/
   #Run CLUMPAK with absolute paths:
   perl CLUMPAK.pl !{customdraw} --id clumpak_run --dir ${base_dir}/output --inputtype admixture --file ${base_dir}/input/!{Qzip} --indtopop ${base_dir}/input/clumpak_pop_map.txt --mclthreshold !{params.mclthreshold} 2> CLUMPAK_!{params.plink_prefix}.stderr > CLUMPAK_!{params.plink_prefix}.stdout
   popd
   #Since there are a ton of different outputs, let's just throw them into a tarball for export:
   tar -czf CLUMPAK_!{params.plink_prefix}_output.tar.gz output/
   '''
}

workflow {
   //Set up file channels for the PLINK inputs:
   plink_bed = file(params.plink_prefix+".bed", checkIfExists: true)
   plink_bim = file(params.plink_prefix+".bim", checkIfExists: true)
   plink_fam = file(params.plink_prefix+".fam", checkIfExists: true)

   //Set up a file channel for the population map:
   pop_map = file(params.pop_map, checkIfExists: true)

   //Set up a file channel for the samples to exclude:
   params.samples_to_exclude = "<(echo)"
   excludesamples = file(params.samples_to_exclude, checkIfExists: true)

   //Set up a file channel for the CLUMPAK install zip:
   CLUMPAK = file(params.clumpak_zip, checkIfExists: true)

   //Set up an optional file channel for distruct's drawing parameters:
   drawparams = file(params.drawparams)

   //Set up value channels for the Ks and replicates:
   minK = params.minK
   maxK = params.maxK
   Channel
      .of(minK..maxK)
      .tap { Ks }
      .subscribe { println "Adding K=${it} to Ks channel" }
   Channel
      .of(1..params.replicates)
      .tap { replicates }
      .subscribe { println "Adding replicate ${it} to replicates channel" }
   //Set up the list of PRNG seeds:
   Random prng = new Random(params.prng_seed)
   Channel
      .fromList(prng.ints(params.replicates, 0, 2147483647).toArray().toList())
      .tap { admixture_seeds }
      .subscribe { println "Adding PRNG seed ${it} to admixture_seeds channel" }

   //Exclude any samples requested from the PLINK files and the pop_map:
   subset(plink_bed, plink_bim, plink_fam, pop_map, excludesamples)

   //Run ADMIXTURE for all Ks and replicates with the subset PLINK files:
   admixture_input = subset.out.plink
      .combine(Ks)
      .combine(replicates.merge(admixture_seeds))
   admixture(admixture_input)

   //Collect all the Qs, rearrange them to match the subset.out.pop_map order, and package them for CLUMPAK:
   Qs = admixture.out.Q
      .collectFile() { [ "admixture_q_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }
   zip_Q(Qs, subset.out.fam, subset.out.pop_map)

   //Run CLUMPAK:
   clumpak(CLUMPAK, subset.out.pop_map, zip_Q.out.zip, drawparams)
}

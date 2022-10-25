#!/usr/bin/env nextflow
/* Pipeline for parallelizing ADMIXPIPE across K and replicates          *
 * Steps:                                                                *
 *  Run ADMIXPIPE per K and per replicate ->                             *
 *  Zip outputs per K, then zip the zips together ->                     *
 *  Run CLUMPAK in a very particular way                                 */

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

//Set up file channels for the PLINK inputs:
plink_bed = file(params.plink_prefix+".bed", checkIfExists: true)
plink_bim = file(params.plink_prefix+".bim", checkIfExists: true)
plink_fam = file(params.plink_prefix+".fam", checkIfExists: true)

//Set up a file channel for the population map:
pop_map = file(params.pop_map, checkIfExists: true)

//Set up a file channel for the CLUMPAK install zip:
CLUMPAK = file(params.clumpak_zip, checkIfExists: true)

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
admixture_seeds = Channel.value(prng.ints(params.replicates, 0, 2147483647).toArray().toList())

//Defaults for cpus, memory, and time for each process:
//ADMIXPIPE ADMIXTURE step:
params.admixture_cpus = 1
params.admixture_mem = 16
params.admixture_timeout = '24h'
//Zipping of Q files for CLUMPAK:
params.zip_q_cpus = 1
params.zip_q_mem = 4
params.zip_q_timeout = '12h'
//CLUMPAK:
params.clumpak_cpus = 1
params.clumpak_mem = 16
params.clumpak_timeout = '24h'

process admixture {
   tag "K${K} r${replicate}"

   cpus params.admixture_cpus
   memory { params.admixture_mem.plus(task.attempt.minus(1).multiply(32))+' GB' }
   time { task.attempt == 2 ? '48h' : params.admixture_timeout }
   errorStrategy { task.exitStatus in ([1]+(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/ADMIXTURE", mode: 'copy', pattern: '*.{P,Q}'

   input:
   val seeds from admixture_seeds
   path plink_bed
   path plink_bim
   path plink_fam
   path pop_map
   each K from Ks
   each replicate from replicates

   output:
   tuple path("admixture_K${K}_r${replicate}.stderr"), path("admixture_K${K}_r${replicate}.stdout") into admixture_logs
   path("${plink_prefix}.${K}_${replicate}.P") into admixture_P
   path("${plink_prefix}.${K}_${replicate}.Q") into admixture_Q
   path("admixture_K${K}_r${replicate}.stdout") into admixture_CV

   shell:
   prng_seed = seeds[replicate.minus(1)]
   plink_prefix = plink_bed.getBaseName()
   '''
   module load !{params.mod_admixture}
   admixture -j!{task.cpus} --seed=!{prng_seed} --cv=!{params.cv_folds} !{plink_prefix}.bed !{K} 2> admixture_K!{K}_r!{replicate}.stderr > admixture_K!{K}_r!{replicate}.stdout
   mv !{plink_prefix}.!{K}.P !{plink_prefix}.!{K}_!{replicate}.P
   mv !{plink_prefix}.!{K}.Q !{plink_prefix}.!{K}_!{replicate}.Q
   '''
/*   module load !{params.mod_miniconda}
   conda activate gaia
   ./admixturePipeline.py -m !{pop_map} -b !{plink_prefix} -k !{K} -K !{K} -n !{task.cpus} -R 1 2> ADMIXPIPE_K!{K}_r!{replicate}.stderr > ADMIXPIPE_K!{K}_r!{replicate}.stdout*/
}

process zip_Q {
   cpus params.zip_q_cpus
   memory { params.zip_q_mem.plus(task.attempt.minus(1).multiply(1))+' GB' }
   time { task.attempt == 2 ? '12h' : params.zip_q_timeout }
   errorStrategy { task.exitStatus in ([1]..(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   input:
   path("admixture_q_paths.tsv") from admixture_Q.collectFile() { [ "admixture_q_paths.tsv", it.getSimpleName()+'\t'+it.getName()+'\t'+it+'\n' ] }

   output:
   path("${params.run_name}.zip") into Q_zip

   shell:
   '''
   while read -a a;
      do
      ln -s ${a[2]} ${a[1]};
   done < admixture_q_paths.tsv
   for K in {!{params.minK}..!{params.maxK}};
      do
      zip !{params.run_name}_K${K}.zip $(find . -name "*.Q" -print | awk -v "K=${K}" '{n=split($1, a, ".");m=split(a[n-1], b, "_");if (b[1] == K) {printf " %s", $1;};}')
   done
   zip !{params.run_name}.zip !{params.run_name}_K*.zip
   '''
}

process clumpak {
   cpus params.clumpak_cpus
   memory { params.clumpak_mem.plus(task.attempt.minus(1).multiply(1))+' GB' }
   time { task.attempt == 2 ? '12h' : params.clumpak_timeout }
   errorStrategy { task.exitStatus in ([1]..(134..140).collect()) ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/ADMIXTURE", mode: 'copy', pattern: '*.tar.gz'

   input:
   path CLUMPAK
   path pop_map
   path Qzip from Q_zip

   output:
   tuple path("${params.clumpak_path}/CLUMPAK_${params.run_name}.stderr"), path("${params.clumpak_path}/CLUMPAK_${params.run_name}.stdout") into clumpak_logs
   path("CLUMPAK_${params.run_name}_output.tar.gz") into clumpak_output

   shell:
   '''
   #Stage the CLUMPAK installation:
   unzip !{CLUMPAK}
   chmod u+x !{params.clumpak_path}/CLUMPP/CLUMPP
   chmod u+x !{params.clumpak_path}/distruct/distruct1.1
   chmod u+x !{params.clumpak_path}/mcl/bin/*
   #Set up the input directory and files:
   mkdir input
   ln -s ../!{pop_map} input/
   ln -s ../!{Qzip} input/
   #Set up the output directory:
   mkdir output
   #Capture the absolute path of the base directory so we can indicate absolute paths in the call:
   base_dir=$(pwd)
   #Set the CWD to the staged CLUMPAK install:
   pushd !{params.clumpak_path}/
   #Run CLUMPAK with absolute paths:
   perl CLUMPAK.pl --id !{params.run_name} --dir ${base_dir}/output --inputtype admixture --file ${base_dir}/input/!{Qzip} --indtopop ${base_dir}/input/!{pop_map} --mclthreshold !{params.mclthreshold} 2> CLUMPAK_!{params.run_name}.stderr > CLUMPAK_!{params.run_name}.stdout
   popd
   #Since there are a ton of different outputs, let's just throw them into a tarball for export:
   tar -czf CLUMPAK_!{params.run_name}_output.tar.gz output/
   '''
}

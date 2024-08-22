# How to run ADMIXTURE.nf on your dataset

## Installation:

`ADMIXTURE.nf` is a Nextflow pipeline written in DSL2, so you will of course need Nextflow installed
in order to run the pipeline. Keep in mind that Nextflow pipelines are simply script files, so the
only installation involved is installation of any necessary dependencies. I haven't written this
pipeline to use containers, so these setup instructions are based around installing dependencies
locally and using the `module` system with [modulefile](https://modules.readthedocs.io/en/stable/modulefile.html)s in order to manage dependency loading.

### Dependencies:

- [nextflow](https://nextflow.io/) (tested with version 22.10.7b5853)
- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) (tested with 1.90 beta 7 and 1.90 beta 6.26)
- [ADMIXTURE 1.3.0](https://dalexander.github.io/admixture/download.html)
- [CLUMPAK](https://tau.evolseq.net/clumpak/download.html) (tested with version from March 26, 2015)
- Perl 5 (tested with 5.32.0)
- Ghostscript (tested with 9.53.3)

PLINK, ADMIXTURE, and CLUMPAK will need modulefiles set up and accessible by the cluster. These
module paths will need to be specified in the `.config` file for the Nextflow run (see the `profile`
section of the example `ADMIXTURE.config`).

CLUMPAK itself has several dependencies pre-packaged with it (e.g. CLUMPP, mcl, distruct), as well
as quite a few Perl module dependencies. Here's a relatively quick way to install these using cpanminus:

```bash
#Load a Perl module if there is no base Perl install or there are issues with it.
#e.g. on Yale's mccleary cluster:
module load Perl/5.32.0-GCCcore-10.2.0

#Install cpanminus, which simplifies Perl module installation and allows for local installs.
#Local installs are usually necessary on shared-usage high-performance computing clusters.
curl -L https://cpanmin.us | perl - App::cpanminus

#Set up the local::lib module to avoid repeated warnings about non-root installations:
cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)

#Install CLUMPAK Perl module dependencies and check that they work properly:
clumpak_deps=("Archive::Zip" "Getopt::Long" "Getopt::Std" "File::Basename" "File::Path" "File::Slurp" "List::Util" "List::MoreUtils" "List::Permutor" "PDF::API2" "PDF::Table" "GD::Graph" "GD::Graph::lines" "GD::Graph::Data" "Scalar::Util" "Statistics::Distributions" "Archive::Extract" "Data::PowerSet" "Array::Utils")
for m in clumpak_deps;
   do
   cpanm ${m};
   perl -M${m} -e 1;
done
```

Keep note if you see any lines containing `Can't locate _____ in @INC`, as that means something in the
Perl module installation failed, or else the `PERL5LIB` environment variable isn't set correctly.

Also, if you get an error about `Can't locate ExtUtils/Manifest.pm in @INC` while installing cpanminus,
there's likely something messed up with your Perl install.

Once these Perl modules are installed, you'll want to set up a modulefile for CLUMPAK that looks
something like this:

```
#%Module1.0#######################################################################
set software "CLUMPAK"
set version "26_03_2015"
#Release 26_03_2015 is latest as of 2023/12/12

proc ModulesHelp { } {
        global version

        puts stderr "\tThis module loads ${software} version ${version}."
}

module-whatis   "Loads ${software} version ${version}"

module load Perl/5.32.0-GCCcore-10.2.0
module load Ghostscript/9.53.3-GCCcore-10.2.0
prepend-path PERL5LIB $::env(HOME)/perl5/lib/perl5/
prepend-path PERL5LIB .
```

You'll have to adjust the Perl and Ghostscript `module load` calls appropriately (or remove if
your base environment already has working versions of both), but keep the `PERL5LIB` lines.
These lines ensure that the Perl modules you installed can be detected, and that the Perl
modules packaged with CLUMPAK are also detected.

## How to run `ADMIXTURE.nf`:

Once you've installed all the dependencies and determined the appropriate modules for them,
you'll need to set up a `.config` file for your run. I've provided an example in `ADMIXTURE.config`
in this repository. There are four general sections to the `.config` file:

1. "Executor" setup (i.e. settings for the cluster)
2. Pipeline parameters (e.g. input globs, output prefixes, thresholds, etc.)
3. Profiles (i.e. parameters specific to a particular cluster or reference)
4. Process settings (e.g. number of cores, memory, timeout, queue to use)

You might also notice a `trace` block before the profiles section. You likely won't need
to adjust this block. It configures Nextflow to save information about jobs to a special
trace file. This information includes whether the process completed, failed, was aborted,
or was not rerun (i.e. "CACHED") when `--resume` was used, as well as the elapsed time,
memory used, exit code, etc. It is a very useful file for doing accounting of the pipeline.

### "Executor" setup:

The "executor" setup section is set up for SLURM, and also specifies the maximum number of
simultaneous queued jobs (`executor.queueSize=300` in the example)
and rate of queueing jobs (`executor.submitRateLimit=100/5min` for 100 every 5 minutes
in the example). This ensures that you can parallelize well, while not spamming the queue
and making your sysadmin's job harder. Ask your cluster sysadmin if you should change these
defaults. I've also added `process.clusterOptions="--requeue"` to the example, since Nextflow
can get thrown off when a SLURM node fails. This option tells SLURM to automatically requeue
these jobs onto working nodes.

### Pipeline parameters:

In Nextflow, variables that can be defined in the config or globally in the pipeline but
are accessible within the processes are typically prefixed with `params.`. In this section
of the `.config` file, we provide variables that determine the location for pipeline outputs,
pipeline inputs, and other pipeline parameters:

Pipeline input variables:

```
params.run_name = "PIBv1"
params.output_prefix = "/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes"
params.output_dir = "${params.output_prefix}/${params.run_name}_ADMIXTURE"
```

This block defines the absolute path where you want all the output directories. In the
example, this would be `/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes/PIBv1_ADMIXTURE`.

```
//Default paths, globs, and regexes:
//Prefix for PLINK input files:
params.plink_prefix = "PIBv1_autosomes_globalMAFge0.01_bSNPs_Choin"
//Path to population map file:
params.pop_map = "admix_pop_map.txt"
```

This block defines the prefix for the PLINK files used as input. The PLINK files should
be in the current working directory whence you execute the pipeline, so they should not
contain any path elements, just a file prefix. The block also defines the path to the
sample-to-population map file, which is a two-column tab-separated file without any header
that has the samples ordered in the order desired for the plotting output. The first
column is the sample ID, and the second column is the population that sample belongs to.


```
//Note that this PRNG seed is actually a kind of "hyperseed", seeding the PRNG that generates the per-replicate seeds for ADMIXTURE:
params.prng_seed = 42
//Range of K to run:
params.minK = 2
params.maxK = 12
//Number of replicates per K:
params.replicates = 20
//Number of folds for cross-validation:
params.cv_folds = 20
//Threshold to use for MCL clustering of components across K values:
params.mclthreshold = 0.9
```

This block contains the general parameters for the pipeline, such as the PRNG hyperseed
to use when generating the seeds for each ADMIXTURE replicate, as well as the range of
K values to test, the number of replicates, the number of folds to use for cross-validation
in ADMIXTURE, and the threshold to use for MCL clustering in CLUMPAK.

```
//Short string to describe the MAF filter being used in this run:
params.MAF_name = "globalMAFge0.01"
//Short string to describe the LD pruning parameters being used in this run:
params.LD_name = "Choin"
```

This block contains variables used in the naming of the output files, specifically
describing the MAF filter threshold and LD pruning approaches used.

### Profiles:

In the example `ADMIXTURE.config` file, the profiles section specifies parameters
specific to the cluster being used, such as module paths, and the path to the CLUMPAK
zip file.

```
profiles {
   mccleary {
      params.mod_plink = "plink/1.90b7"
      params.mod_admixture = "ADMIXTURE/1.3.0"
      params.mod_clumpak = "CLUMPAK/26_03_2015"
      //Path to CLUMPAK install zip:
      params.clumpak_zip = "/home/pfr8/bin/CLUMPAK/20150326/26_03_2015_CLUMPAK.zip"
      //Relative path to CLUMPAK Perl modules:
      params.clumpak_path = "26_03_2015_CLUMPAK/CLUMPAK"
      process.queue = 'day'
   }
}
```

It also contains the relative path to the CLUMPAK scripts and Perl modules, as well
as a specification of the default queue to use on the mccleary cluster. On previous
clusters, this may have been something like `process.queue = 'general'` instead.

### Process settings:

This section is fairly repetitive, but specifies the SLURM job submission parameters
for each Nextflow process.

```
//SLURM submission parameters:
//Memory specified in GiB unless otherwise indicated
//Defaults for cpus, memory, and time for each process:
//Input subsetting:
//Max 1 retry
//Adds 16 GiB per retry
//Timeout to 24h on retry
params.subset_cpus = 4
params.subset_mem = 3
params.subset_timeout = '1h'
params.subset_queue = 'day'
```

So this just specifies that the input subsetting process should be allocated 4 cores,
3 GiB of memory, and 1 hour timeout, and submitted on the 'day' partition/queue.

```
//ADMIXTURE step:
//Max 1 retry
//Adds 32 GiB per retry
//Timeout to params.admixture_long_timeout on retry
params.admixture_cpus = 4
params.admixture_mem = 16
params.admixture_timeout = '72h'
params.admixture_long_timeout = '168h'
params.admixture_queue = 'week'
params.admixture_long_queue = 'week'
```

In the case of ADMIXTURE, it's possible that the job might fail and get retried with
a longer timeout. Hence, I distinguish between `params.admixture_queue` and
`params.admixture_long_queue`, the latter of which is used as the partition/queue
upon Nextflow retrying submission. It looks a bit weird here, since both the initial
queue and the long queue are 'week', but if `params.admixture_timeout = '24h'`, then
I would set `params.admixture_queue = 'day'`, and you'd see why it's necessary.

### Executing the Nextflow pipeline:

Once you have your `.config` file set up, make sure you have the necessary inputs
in place and where the `.config` file says they should be. These inputs are:

1. The input PLINK files (`.bed`, `.bim`, and `.fam`)
2. The sample-to-population map TSV
3. The optional `params.samples_to_exclude` file of sample IDs to exclude

Once those files are in place, be sure to create a scratch directory for the
pipeline run somewhere where you have a decent amount of free space. Nextflow
will create a bunch of directories in the scratch directory in a hierarchy
where the first level of directories are the two hexadecimal digit prefixes
of the Nextflow jobs, and within each first level directory is one directory
per job that has that prefix.

Once you have all of this in place, you can execute the pipeline (with Nextflow loaded) using:
```bash
NXF_OPTS="-Xms500M -Xmx900M" /usr/bin/time -v nextflow -c [path to .config file] -bg run [path to .nf file] -profile [profile] -w [path to scratch directory] 2> [STDERR file for this run] > [STDOUT file for this run]
```

You can watch the pipeline run by then running:
```bash
tail -n+1 -f [STDOUT of the Nextflow command] | less
```

This looks a bit different from a more basic Nextflow command for a few reasons. First,
the `NXF_OPTS` environment variable provides the Java options controlling memory usage
of the Nextflow process running on the head node of the cluster. It sets the minimum
memory allocation to 500 MB, and the upper limit at 900 MB. I set this because Yale
clusters generally automatically kill any processes on the head node using more than
1 GB of memory. Nextflow pipelines generally don't use a huge amount of memory unless
there are a *ton* of simultaneous jobs or an incredibly complex DAG of jobs. Second,
the `-bg` option makes Nextflow run in the background, so it won't be killed by a SIGPIPE
if your SSH connection disconnects while the pipeline is running. You can safely logout
and the pipeline will continue running. `-profile [profile]` tells Nextflow to use the
parameters specified in the `.config` file listed under the `[profile]` profile, so
`-profile mccleary` would be what we want for the example `ADMIXTURE.config` file.
`-w` specifies the scratch directory, and `-c` specifies the `.config` file. The rest
is simply redirecting STDERR and STDOUT so that you can safely log out and maintain
a log of the run.

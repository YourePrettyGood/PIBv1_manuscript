# Processing and analysis pipelines for the manuscript titled "Long-term isolation and archaic introgression shape functional genetic variation in Near Oceania"

## Repository organization

This repository is organized into Processing Pipelines, Analysis Pipelines, Plotting, Configs, and Instructions subdirectories.

### [Processing pipelines](/PROCESSING.md)

The processing pipelines were used for initial QC and screening of
newly-sequenced Near Oceanic samples, mapping and variant calling of new and
publicly available human whole genomes, joint calling of the full dataset,
and phasing of the genotypes. These pipelines are provided as Nextflow DSL1
pipelines (they were developed before DSL2 became default) and are authored
by Patrick F. Reilly.

### [Analysis pipelines](/ANALYSIS.md)

The analysis pipelines were used for any analysis of the jointly called and
filtered dataset. Many analysis pipelines are provided as Nextflow pipelines
(mostly DSL1, some as DSL2), while some analyses performed by other co-authors
are provided as submodules pointing to their respective repositories.
The [`ANALYSIS.md`](/ANALYSIS.md) document provides a lot of details about
the inputs, commands, thought process, and choices involved in each analysis,
so please read that for details.

Analyses include:

- Joint callset filtering and filter selection (DSL1: `filter_VCF.nf`)
- Calculation of VCF statistics (including mutation matrix, Ts/Tv, per-sample stats, Mendelian error rate from trios, genotype concordance with array data) (DSL1: `VCF_stats.nf`)
- Preparation of MAF-filtered and LD-pruned PLINK files for other analyses (DSL1: `MAFLD_VCF.nf`)
- Preliminary relatedness and contamination analyses (DSL1: `BAM_contam_QC.nf`)
- Novel and population-private variant analysis
- ADMIXTURE analysis (DSL1: `ADMIXPIPE.nf`, DSL2: `ADMIXTURE.nf`)
- Principal components analysis (PCA as an R script: `SNPRelate_analyses.R`)
- Population structure analyses with ADMIXTOOLS (outgroup-f<sub>3</sub>, f<sub>4</sub>)
- Local ancestry inference with RFMix (see Chang's repository)
- ROH analysis (see Audrey's repository)
- Demographic inference (SMC++ in DSL1: `smcpp.nf`, MSMC2 including cross-coalescence and MSMC-IM in DSL2: `msmc2.nf`)
- Archaic introgression (Sprime in DSL1: `sprime.nf`, Core haplotypes in DSL1: `adaptiveintrogression.nf`, phased projections of Sprime tracts in DSL1: `sprime_phased_projections.nf`, ArchaicSeeker2.0 and MultiWaver3.1 in DSL1: `archaicseeker.nf`, and see Chang's repository for D and f<sub>4</sub> ratio statistics as well as PCA projection analyses)
- Gaussian and Beta mixture model analyses of Sprime tract match rates (R script: `SprimeGMMBMM.R` and related scripts for supplementary figures)
- Adaptive introgression scans (see Dani's repository, and the performance evaluation section including the simulation pipeline in DSL2: `adaptive_introgression_simulations.nf`)
- Variant effect prediction by Ensembl VEP (DSL2: `VEP.nf`)
- Functional genomics analyses (see Stephen's repository)
- MPRA analyses (see Stephen's repository)
- Generation of haplotype plots (DSL2: `haplotype_plot.nf`)

Note that the ArchaicSeeker2.0 and MSMC2/MSMC-IM analyses are not included in the current manuscript.
Also, the `callable_masks.nf` pipeline was run with the intent of using the
outputs for ARGWeaver, but that didn't make it into the manuscript.

### [Plotting scripts](/PLOTS.md)

Scripts used for generating most of the main text and supplementary figures
are linked in the [`PLOTS.md`](/PLOTS.md) document. Some figures were
composited and/or reformatted for publication using Adobe Illustrator, so
I've tried to annotate those in the document or in their own scripts.

### Configs

Each Nextflow pipeline relies on a `.config` file to provide run-specific
inputs and thresholds. In the `Configs` directory, we provide both the
original `.config` files used for the manuscript runs (subdirectory `Used`)
as well as example `.config` files (subdirectory `Example`) for broader
adaptation and use.

Modulefiles for dependencies of each pipeline can be found in the `modules`
subdirectory.

### Instructions

This is a work in progress, but each Nextflow pipeline will come along with
documentation in the form of a markdown document in the `Instructions`
subdirectory.

## General usage and dependency information

In general, these Nextflow pipelines will all work with Nextflow 22.10.7,
though make sure to set `NXF_DEFAULT_DSL=1` for any DSL1 pipelines. The
DSL2 pipelines were run with Nextflow 22.10.7b5853 or 24.10.2b9532, while
the DSL1 pipelines were run with either 21.10.5b5658 or 21.10.3b5655.

Dependencies of each pipeline are generally specified by the `params.mod_*`
variables in the `.config` file, though these are set up as modulefiles
rather than Docker/Singularity/AppTainer containers. Exceptions to this
include conda environments (where `params.mod_miniconda` must be specified)
and R packages (where `params.mod_R` must be specified). I try to document
all relevant dependencies in the `Instructions` for each pipeline.
Modulefiles are also provided in the `Configs` directory under the `modules`
subdirectory.

Most pipelines were executed with the general command template:

If run on the head node:

```bash
NXF_OPTS="-Xms500M -Xmx900M" /usr/bin/time -v nextflow -c [path to .config file] run -bg [absolute path to .nf pipeline script] -profile [comma-separated list of profiles] -w [scratch/working directory] 2> [STDERR log] > [STDOUT log]
```

or if run in it's own submitted job:

```bash
NXF_OPTS="-Xms500M -Xmx3500M" /usr/bin/time -v nextflow -c [path to .config file] run -ansi-log false [absolute path to .nf pipeline script] -profile [comma-separated list of profiles] -w [scratch/working directory] 2> [STDERR log] > [STDOUT log]
```

Note that you will need to add `NXF_DEFAULT_DSL=1` as a prefix for any DSL1
pipelines when using Nextflow 22.x.

You can set `-Xmx` to a larger value than `900M` like `3500M`, though some pipelines
shouldn't need more. Increasing `-Xmx` is generally only needed when the
pipeline generates a lot of jobs (e.g. `adaptiveintrogression.nf` or `VEP.nf`).
I used to set it to `900M` because our cluster would OOM kill the nextflow process
at 1 GB of RSS, and an OOM kill of the main nextflow process causes a dirty
pipeline exit that does not resume correctly and truncates the log and cache.

## Limitations

These pipelines were generally written with human data mapped to the hs37d5
reference genome in mind. I tried to write them in a fairly generalized way
so that they're compatible with other human references, and possibly even
non-human species, but cannot guarantee compatibility.

Only a select few pipelines have been adapted and tested for compatibility
with GRCh38:

- `ADMIXTURE.nf`
- `smcpp.nf`
- `sprime.nf`
- `adaptiveintrogression.nf`

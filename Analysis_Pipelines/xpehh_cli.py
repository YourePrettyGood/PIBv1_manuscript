#!/usr/bin/env python3
import argparse
import sys
import time
import allel
import numpy
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description='Command-line interface (CLI) for evaluating XP-EHH, XP-nSL, and iHS with scikit-allel')
    parser.add_argument('-i', '--input_vcf', help='VCF of variants for XP-EHH evaluation', required=True)
    parser.add_argument('-m', '--sample_metadata', help='TSV of sample metadata including columns for sample ID, population, and superpopulation (and a header line with column names)', required=True)
    parser.add_argument('-s', '--sample_id_col', help='Name of column containing sample IDs in the sample metadata file (default: %(default)s)', default='SampleID')
    parser.add_argument('-p', '--pop_col', help='Name of the column containing population labels in the sample metadata file (default: %(default)s)', default='Population')
    parser.add_argument('--superpop_col', help='Name of the column containing superpopulation labels in the sample metadata file')
    parser.add_argument('-g', '--genetic_map', help='Genetic map file (for a single chromosome) in TSV format with a header line and minimally columns for marker physical position ("pos") and marker genetic map position ("cM")', required=True)
    parser.add_argument('--pos_col_index', help='0-based index of the physical position column in the genetic map file', type=int, required=True)
    parser.add_argument('--cm_col_index', help='0-based index of the genetic map position column in the genetic map file', type=int, required=True)
    parser.add_argument('-t', '--target_pop', help='Target population for PBS', required=True)
    parser.add_argument('-r', '--reference_pop', help='Ingroup reference population for PBS', required=True)
    parser.add_argument('-o', '--output_tsv', help='Name/path of TSV to output containing PBS results', required=True)
    return parser.parse_args()

def evaluate_stats(target_haps, ref_haps, pos, gmap, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True):
    """
    """

    #The code below largely comes from Dani's run_xpehh.py, just extended to
    # other stats calculated by scikit-allel from 1 or 2 populations:
    xpehh_raw = allel.xpehh(target_haps, ref_haps, pos, map_pos=gmap, min_ehh=min_ehh, include_edges=include_edges, gap_scale=gap_scale, max_gap=max_gap, is_accessible=is_accessible, use_threads=use_threads)
    xpehh_raw[~numpy.isfinite(xpehh_raw)] = numpy.nan
    xpehh = allel.standardize(xpehh_raw)

    #XP-nSL:
    xpnsl_raw = allel.xpnsl(target_haps, ref_haps, use_threads=use_threads)
    xpnsl_raw[~numpy.isfinite(xpnsl_raw)] = numpy.nan
    xpnsl = allel.standardize(xpnsl_raw)

    #iHS:
    #iHS is a bit funky, as the scikit-allel code doesn't handle
    # multiallelic sites properly, so we have to filter them and reintroduce
    # as NaNs in the output.
    ihs_raw = numpy.empty(target_haps.n_variants)
    ihs_raw.fill(numpy.nan)
    biallelic_sites = target_haps.count_called(axis=1) == target_haps.count_ref(axis=1) + target_haps.count_call(1, axis=1)
    ihs_raw[biallelic_sites] = allel.ihs(target_haps[biallelic_sites], pos[biallelic_sites], map_pos=gmap[biallelic_sites], min_ehh=min_ehh, include_edges=include_edges, gap_scale=gap_scale, max_gap=max_gap, is_accessible=is_accessible, use_threads=use_threads)
    ihs_raw[~numpy.isfinite(ihs_raw)] = numpy.nan
    ihs, ihs_aac_bins = allel.standardize_by_allele_count(ihs_raw, target_haps.count_alleles()[:,1])

    return numpy.column_stack((xpehh, xpnsl, ihs))

def get_pop_lists(metadata_path, id_col, pop_col, superpop_col, target, ref):
    id_pop_map = list()
    with open(metadata_path, 'r') as metadata:
        for i,line in enumerate(metadata):
            cols = line.rstrip().split('\t')
            if i == 0:
                try:
                    id_col_index = cols.index(id_col)
                    pop_col_index = cols.index(pop_col)
                except ValueError as v:
                    print(f"id_col ({id_col}) or pop_col ({pop_col}) missing from metadata file header ({line}), cannot proceed.", file=sys.stderr, flush=True)
                    raise
                if superpop_col is not None:
                    try:
                        superpop_col_index = cols.index(superpop_col)
                    except ValueError as v:
                        print(f"superpop_col ({superpop_col}) missing from metadata file header ({line}), cannot proceed.", file=sys.stderr, flush=True)
                        raise
            else:
                if superpop_col is not None:
                    id_pop_map.append((cols[id_col_index], cols[pop_col_index], cols[superpop_col_index]))
                else:
                    id_pop_map.append((cols[id_col_index], cols[pop_col_index]))

    target_pops = target.split(',')
    ref_pops = ref.split(',')

    if superpop_col is not None:
        target_ids = [id for id, pop, superpop in id_pop_map if pop in target_pops or superpop in target_pops]
        ref_ids = [id for id, pop, superpop in id_pop_map if pop in ref_pops or superpop in ref_pops]
    else:
        target_ids = [id for id, pop in id_pop_map if pop in target_pops]
        ref_ids = [id for id, pop in id_pop_map if pop in ref_pops]

    if len(target_ids) == 0 or len(ref_ids) == 0:
        raise ValueError(f'No samples found in metadata for at least one of the groups required for XP-EHH. Target group was {target} and reference group was {ref}.')

    return (target_ids, ref_ids)

def interpolate_gmap(gmap_path, pos, pos_col_index, cM_col_index):
    gmap_phys, gmap_gen = numpy.loadtxt(gmap_path, delimiter='\t', skiprows=1, usecols=(pos_col_index, cM_col_index), unpack=True)
    if numpy.isclose(gmap_phys, gmap_phys[0]).all() or numpy.isclose(gmap_gen, gmap_gen[0]).all():
        raise ValueError(f'Physical or genetic map positions of all elements in the genetic map ({gmap_path}) are too similar to each other. Did you misspecify --pos_col_index ({pos_col_index}) or --cm_col_index ({cM_col_index}) for this genetic map file?')
    return numpy.interp(pos, gmap_phys, gmap_gen)

def write_xpehh(out_path, sites, xpehh_arr, target, ref):
    #Double-check that sites and xpehh_arr are of the same length:
    if len(sites) != xpehh_arr.shape[0]:
        raise ValueError(f'Somehow the sites list has different length from the XP-EHH statistic array: {len(sites)} != {xpehh_arr.shape[0]}')
    suffix = f'{target}\t{ref}'
    with (gzip.open if out_path.endswith('.gz') else open)(out_path, 'wt') as xpehh_out:
        #Print a header line describing the columns:
        print('#CHROM\tPOS\tCM\tID\tXPEHH\tXPNSL\tIHS\tTARGET\tREF', file=xpehh_out)
        #Now iterate through the windows:
        for i in range(xpehh_arr.shape[0]):
            chrom, pos, cm, id = sites[i]
            xpehh, xpnsl, ihs = xpehh_arr[i,]
            print(f'{chrom}\t{pos:d}\t{cm:g}\t{id}\t{xpehh:g}\t{xpnsl:g}\t{ihs:g}\t{suffix}', file=xpehh_out)

def run_xpehh(args):
    #Determine start time and CPU time:
    t_0 = time.perf_counter_ns()
    tcpu_0 = time.process_time_ns()

    #Parse the sample ID lists out from the sample metadata file:
    target, ref = get_pop_lists(args.sample_metadata, args.sample_id_col, args.pop_col, args.superpop_col, args.target_pop, args.reference_pop)

    t_poplist = time.perf_counter_ns()
    tcpu_poplist = time.process_time_ns()
    elapsed_poplist = (t_poplist - t_0) / 1_000_000_000
    elapsedcpu_poplist = (tcpu_poplist - tcpu_0) / 1_000_000_000

    print(f'Extracted relevant sample lists from sample metadata file {args.sample_metadata} in {elapsed_poplist} seconds ({elapsedcpu_poplist} CPU-seconds)', file=sys.stderr, flush=True)

    #Generate the overall sample ID list for extracting from the VCF:
    samples = set(target).union(ref)
    if len(samples) < len(target) + len(ref):
        raise ValueError(f'There appears to be some overlap in samples between target group {target} and reference group {ref}. These should be mutually exclusive groups.')

    #Read in the VCF, but only the relevant parts:
    vcf_fields = ["CHROM", "POS", "ID", "GT", "samples"]
    vcf = allel.read_vcf(args.input_vcf, samples=samples, fields=vcf_fields)

    #Before we get too far, check that the VCF only contains one chromosome,
    # since I haven't implemented interpolate_gmap() for multiple chromosomes
    # yet:
    if numpy.unique(vcf['variants/CHROM']).shape[0] > 1:
        raise ValueError(f'The input VCF {args.input_vcf} appears to contain variants from multiple chromosomes. Evaluating XP-EHH and related statistics on multi-chromosome VCFs is not currently implemented.')

    t_readvcf = time.perf_counter_ns()
    tcpu_readvcf = time.process_time_ns()
    elapsed_readvcf = (t_readvcf - t_poplist) / 1_000_000_000
    elapsedcpu_readvcf = (tcpu_readvcf - tcpu_poplist) / 1_000_000_000

    print(f'Parsed VCF {args.input_vcf} in {elapsed_readvcf} seconds ({elapsedcpu_readvcf} CPU-seconds)', file=sys.stderr, flush=True)

    #Extract the genotype array so we can subset out the three groups and count alleles:
    gts = allel.GenotypeArray(vcf['calldata/GT'])

    t_gtarr = time.perf_counter_ns()
    tcpu_gtarr = time.process_time_ns()
    elapsed_gtarr = (t_gtarr - t_readvcf) / 1_000_000_000
    elapsedcpu_gtarr = (tcpu_gtarr - tcpu_readvcf) / 1_000_000_000

    print(f'Generated allel.GenotypeArrays in {elapsed_gtarr} seconds ({elapsedcpu_gtarr} CPU-seconds)', file=sys.stderr, flush=True)

    #Generate the index lists from the VCF object's sample list:
    target_samples = [i for i, id in enumerate(vcf['samples']) if id in target]
    ref_samples = [i for i, id in enumerate(vcf['samples']) if id in ref]

    #Do the subsetting and extraction of phased haplotypes
    target_haps = allel.HaplotypeArray(gts.subset(sel1=target_samples).haploidify_samples(), dtype='i1')
    ref_haps = allel.HaplotypeArray(gts.subset(sel1=ref_samples).haploidify_samples(), dtype='i1')

    t_haps = time.perf_counter_ns()
    tcpu_haps = time.process_time_ns()
    elapsed_haps = (t_haps - t_gtarr) / 1_000_000_000
    elapsedcpu_haps = (tcpu_haps - tcpu_gtarr) / 1_000_000_000

    print(f'Extracted haplotypes in {elapsed_haps} seconds ({elapsedcpu_haps} CPU-seconds)', file=sys.stderr, flush=True)

    #Extract the sites to evaluate XP-EHH and other stats and interpolate
    # genetic map positions:
    pos = vcf['variants/POS']
    gmap = interpolate_gmap(args.genetic_map, pos, args.pos_col_index, args.cm_col_index)
    sites = list(zip(vcf['variants/CHROM'], pos, gmap, vcf['variants/ID']))

    t_sites = time.perf_counter_ns()
    tcpu_sites = time.process_time_ns()
    elapsed_sites = (t_sites - t_haps) / 1_000_000_000
    elapsedcpu_sites = (tcpu_sites - tcpu_haps) / 1_000_000_000

    print(f'Determined sites for XP-EHH scan in {elapsed_sites} seconds ({elapsedcpu_sites} CPU-seconds)', file=sys.stderr, flush=True)

    #Evaluate XP-EHH and other statistics:
    xpehh = evaluate_stats(target_haps, ref_haps, pos, gmap)

    t_xpehh = time.perf_counter_ns()
    tcpu_xpehh = time.process_time_ns()
    elapsed_xpehh = (t_xpehh - t_sites) / 1_000_000_000
    elapsedcpu_xpehh = (tcpu_xpehh - tcpu_sites) / 1_000_000_000

    print(f'Calculated XP-EHH and other statistics for all sites in {elapsed_xpehh} seconds ({elapsedcpu_xpehh} CPU-seconds)', file=sys.stderr, flush=True)

    #Generate the output file:
    write_xpehh(args.output_tsv, sites, xpehh, args.target_pop, args.reference_pop)

    t_writexpehh = time.perf_counter_ns()
    tcpu_writexpehh = time.process_time_ns()
    elapsed_writexpehh = (t_writexpehh - t_xpehh) / 1_000_000_000
    elapsedcpu_writexpehh = (tcpu_writexpehh - tcpu_xpehh) / 1_000_000_000
    elapsed_total = (t_writexpehh - t_0) / 1_000_000_000
    elapsedcpu_total = (tcpu_writexpehh - tcpu_0) / 1_000_000_000

    print(f'Output XP-EHH and other stats for each site to {args.output_tsv} in {elapsed_writexpehh} seconds ({elapsedcpu_writexpehh} CPU-seconds)', file=sys.stderr, flush=True)
    print(f'Total XP-EHH run took {elapsed_total} seconds ({elapsedcpu_total} CPU-seconds)', file=sys.stderr, flush=True)

if __name__=="__main__":
    args = parse_args()
    run_xpehh(args)

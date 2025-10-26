#!/usr/bin/env python3
import argparse
import sys
import time
import allel
import numpy
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description='Command-line interface (CLI) for evaluating PBS-related statistics with scikit-allel')
    parser.add_argument('-i', '--input_vcf', help='VCF of variants for PBS evaluation', required=True)
    parser.add_argument('-m', '--sample_metadata', help='TSV of sample metadata including columns for sample ID and population (and a header line with column names)', required=True)
    parser.add_argument('-s', '--sample_id_col', help='Name of column containing sample IDs in the sample metadata file (default: %(default)s)', default='SampleID')
    parser.add_argument('-p', '--pop_col', help='Name of the column containing population labels in the sample metadata file (default: %(default)s)', default='Population')
    parser.add_argument('--superpop_col', help='Name of the column containing superpopulation labels in the sample metadata file')
    parser.add_argument('-t', '--target_pop', help='Target population for PBS', required=True)
    parser.add_argument('-r', '--reference_pop', help='Ingroup reference population for PBS', required=True)
    parser.add_argument('-e', '--outgroup_pop', help='Outgroup/exclusionary population for PBS', required=True)
    parser.add_argument('-o', '--output_tsv', help='Name/path of TSV to output containing PBS results', required=True)
    parser.add_argument('-w', '--window_size', type=int, help='Size of the window (in number of variants) for calculating PBS (default: %(default)d)', default=1)
    parser.add_argument('--window_start', type=int, help='0-based index (i.e. NOT position) of the first variant to consider for windowing (default: %(default)d)', default=0)
    parser.add_argument('--window_end', type=int, help='0-based index (i.e. NOT position) of the last variant to consider for windowing (default: the last variant available)')
    parser.add_argument('--window_step', type=int, help='Window step size (in number of variants) for calculating PBS (default: same size as the window)')
    parser.add_argument('--no_pbsn1', action='store_true', help='Calculate the original PBS statistic of Yi et al. 2010 rather than the normalized PBSn1 statistic of Malaspinas et al. 2016 (default: calculate PBSn1)')
    parser.add_argument('--clip_neg_pbs', action='store_true', help='Set any negative PBSn1/PBS scores to 0')
    return parser.parse_args()

def evaluate_stats(target_AC, ref_AC, og_AC, window_size, window_start=0, window_stop=None, window_step=None, clip_neg_pbs=None):
    """
    """

    #The code below largely comes from scikit-allel's allel.pbs() function source, but is adapted to output
    # several related statistics simultaneously, like PBS, PBSn1, LSBL, Akey's d_i, etc.
    #Cast allele counts and check for compatibility as done in scikit-allel:
    ac1 = allel.model.ndarray.AlleleCountsArray(target_AC)
    ac2 = allel.model.ndarray.AlleleCountsArray(ref_AC)
    ac3 = allel.model.ndarray.AlleleCountsArray(og_AC)
    allel.util.check_dim0_aligned(ac1, ac2, ac3)

    #Compute the three pairwise F_{ST} vectors using the Hudson, Slatkin, and Maddison 1992 estimator:
    fst12 = allel.stats.fst.moving_hudson_fst(ac1, ac2, size=window_size, start=window_start, stop=window_stop, step=window_step)
    fst13 = allel.stats.fst.moving_hudson_fst(ac1, ac3, size=window_size, start=window_start, stop=window_stop, step=window_step)
    fst23 = allel.stats.fst.moving_hudson_fst(ac2, ac3, size=window_size, start=window_start, stop=window_stop, step=window_step)

    #Clip F_{ST} values to avoid F_{ST}==1.0, since that would result in infinities when calculating the T vectors for PBS:
    for x in fst12, fst13, fst23:
        numpy.clip(x, a_min=0, a_max=0.99999, out=x)

    #Calculate the three tree distance estimates for PBS:
    t12 = -numpy.log(1 - fst12)
    t13 = -numpy.log(1 - fst13)
    t23 = -numpy.log(1 - fst23)

    #Calculate the PBS estimate of Yi et al. 2010:
    pbs = (t12 + t13 - t23) / 2

    #Normalize PBS into PBSn1 from Malaspinas et al. 2016 Nature:
    pbsn1_denom = 1 + (t12 + t13 + t23) / 2
    pbsn1 = pbs / pbsn1_denom

    #Calculate the LSBL estimate of Shriver et al. 2004 Human Genomics:
    lsbl = (fst12 + fst13 - fst23) / 2

    #Calculate PBE from Yassin et al. 2016 PNAS:
    pbe = pbs - t23 * numpy.nanmedian(pbs) / numpy.nanmedian(t23)

    #Calculate d_i from Akey et al. 2010 PNAS:
    #Slight liberties taken in interpretation here, as we're only handling 3 populations here and only
    # evaluating d_i with respect to the first population, so the sum over j != i just means we include
    # fst12 and fst13 but not fst23. It might be beneficial to include more populations in the sum, but
    # only three are included here for simplicity of comparison with the other PBS-like stats.
    di = (fst12 - numpy.nanmean(fst12)) / numpy.nanstd(fst12) + (fst13 - numpy.nanmean(fst13)) / numpy.nanstd(fst13)

    #Clip the PBS and PBSn1 values if requested:
    if clip_neg_pbs:
        return numpy.column_stack((fst12, numpy.clip(pbs, 0, None), numpy.clip(pbsn1, 0, None), lsbl, pbe, di))
    else:
        return numpy.column_stack((fst12, pbs, pbsn1, lsbl, pbe, di))

def get_pop_lists(metadata_path, id_col, pop_col, superpop_col, target, ref, og):
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
    og_pops = og.split(',')

    if superpop_col is not None:
        target_ids = [id for id, pop, superpop in id_pop_map if pop in target_pops or superpop in target_pops]
        ref_ids = [id for id, pop, superpop in id_pop_map if pop in ref_pops or superpop in ref_pops]
        og_ids = [id for id, pop, superpop in id_pop_map if pop in og_pops or superpop in og_pops]
    else:
        target_ids = [id for id, pop in id_pop_map if pop in target_pops]
        ref_ids = [id for id, pop in id_pop_map if pop in ref_pops]
        og_ids = [id for id, pop in id_pop_map if pop in og_pops]

    if len(target_ids) == 0 or len(ref_ids) == 0 or len(og_ids) == 0:
        raise ValueError(f'No samples found in metadata for at least one of the groups required for PBS. Target group was {target}, reference group was {ref}, and outgroup/exclusionary group was {og}.')

    return (target_ids, ref_ids, og_ids)

def write_pbs(out_path, windows, pbs_arr, target, ref, og, window_size):
    #Double-check that windows and pbs are of the same length:
    if len(windows) != pbs_arr.shape[0]:
        raise ValueError(f'Somehow the windows list has different length from the PBS statistic array: {len(windows)} != {pbs_arr.shape[0]}')
    suffix = f'{target}\t{ref}\t{og}'
    with (gzip.open if out_path.endswith('.gz') else open)(out_path, 'wt') as pbs_out:
        #Print a header line describing the columns:
        if window_size > 1:
            print('#CHROM\tSTARTPOS\tENDPOS\tFST\tPBS\tPBSn1\tLSBL\tPBE\tDi\tTARGET\tREF\tOG', file=pbs_out)
            #Now iterate through the windows:
            for i in range(pbs_arr.shape[0]):
                chrom, start, end = windows[i]
                fst, pbs, pbsn1, lsbl, pbe, di = pbs_arr[i,]
                print(f'{chrom}\t{start:d}\t{end:d}\t{fst:g}\t{pbs:g}\t{pbsn1:g}\t{lsbl:g}\t{pbe:g}\t{di:g}\t{suffix}', file=pbs_out)
        else:
            print('#CHROM\tPOS\tID\tFST\tPBS\tPBSn1\tLSBL\tPBE\tDi\tTARGET\tREF\tOG', file=pbs_out)
            #Now iterate through the windows:
            for i in range(pbs_arr.shape[0]):
                chrom, pos, id = windows[i]
                fst, pbs, pbsn1, lsbl, pbe, di = pbs_arr[i,]
                print(f'{chrom}\t{pos:d}\t{id}\t{fst:g}\t{pbs:g}\t{pbsn1:g}\t{lsbl:g}\t{pbe:g}\t{di:g}\t{suffix}', file=pbs_out)

def run_pbs(args):
    #Determine start time and CPU time:
    t_0 = time.perf_counter_ns()
    tcpu_0 = time.process_time_ns()

    #Parse the sample ID lists out from the sample metadata file:
    target, ref, og = get_pop_lists(args.sample_metadata, args.sample_id_col, args.pop_col, args.superpop_col, args.target_pop, args.reference_pop, args.outgroup_pop)

    t_poplist = time.perf_counter_ns()
    tcpu_poplist = time.process_time_ns()
    elapsed_poplist = (t_poplist - t_0) / 1_000_000_000
    elapsedcpu_poplist = (tcpu_poplist - tcpu_0) / 1_000_000_000

    print(f'Extracted relevant sample lists from sample metadata file {args.sample_metadata} in {elapsed_poplist} seconds ({elapsedcpu_poplist} CPU-seconds)', file=sys.stderr, flush=True)

    #Generate the overall sample ID list for extracting from the VCF:
    samples = set(target).union(ref, og)
    if len(samples) < len(target) + len(ref) + len(og):
        raise ValueError(f'There appears to be some overlap in samples between target group {target}, reference group {ref}, and outgroup/exclusionary group {og}. These should be mutually exclusive groups.')

    #Read in the VCF, but only the relevant parts:
    vcf_fields = ["CHROM", "POS", "ID", "GT", "samples"]
    vcf = allel.read_vcf(args.input_vcf, samples=samples, fields=vcf_fields)

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
    og_samples = [i for i, id in enumerate(vcf['samples']) if id in og]

    #Do the subsetting and allele counting:
    target_AC = gts.subset(sel1=target_samples).count_alleles()
    ref_AC = gts.subset(sel1=ref_samples).count_alleles()
    og_AC = gts.subset(sel1=og_samples).count_alleles()

    t_ACs = time.perf_counter_ns()
    tcpu_ACs = time.process_time_ns()
    elapsed_ACs = (t_ACs - t_gtarr) / 1_000_000_000
    elapsedcpu_ACs = (tcpu_ACs - tcpu_gtarr) / 1_000_000_000

    print(f'Counted alleles in {elapsed_ACs} seconds ({elapsedcpu_ACs} CPU-seconds)', file=sys.stderr, flush=True)

    if args.window_size > 1:
        #Determine the PBS window characteristics:
        window_indices = allel.index_windows(vcf['variants/POS'], size=args.window_size, start=args.window_start, stop=args.window_end, step=args.window_step)
        #Keep in mind that the interval (window_start, window_stop) yielded by allel.index_windows() is right-exclusive,
        # so we need to subtract 1 from the end_index/window_stop (i.e. it's not actually the window_stop/end_index)...
        windows = [(vcf['variants/CHROM'][start_index], vcf['variants/POS'][start_index], vcf['variants/POS'][end_index-1]) for start_index, end_index in window_indices]
    else:
        #Windows are size 1, so just evaluate PBS on each variant:
        if args.window_start is None:
            start = 0
        else:
            start = max(0, args.window_start)
        if args.window_end is None:
            end = len(vcf['variants/POS'])
        else:
            end = min(args.window_end, len(vcf['variants/POS']))
        windows = list(zip(vcf['variants/CHROM'], vcf['variants/POS'], vcf['variants/ID']))[start:end+1]

    t_window = time.perf_counter_ns()
    tcpu_window = time.process_time_ns()
    elapsed_window = (t_window - t_ACs) / 1_000_000_000
    elapsedcpu_window = (tcpu_window - tcpu_ACs) / 1_000_000_000

    print(f'Determined windows for PBS scan in {elapsed_window} seconds ({elapsedcpu_window} CPU-seconds)', file=sys.stderr, flush=True)

#    #Run PBS:
#    pbs = allel.pbs(target_AC, ref_AC, og_AC, window_size=args.window_size, window_start=args.window_start, window_stop=args.window_end, window_step=args.window_step, normed=not args.no_pbsn1)
    #Evaluate the PBS-related statistics:
    pbs = evaluate_stats(target_AC, ref_AC, og_AC, window_size=args.window_size, window_start=args.window_start, window_stop=args.window_end, window_step=args.window_step, clip_neg_pbs=args.clip_neg_pbs)

    t_pbs = time.perf_counter_ns()
    tcpu_pbs = time.process_time_ns()
    elapsed_pbs = (t_pbs - t_window) / 1_000_000_000
    elapsedcpu_pbs = (tcpu_pbs - tcpu_window) / 1_000_000_000

    print(f'Calculated PBS and related statistics for all windows in {elapsed_pbs} seconds ({elapsedcpu_pbs} CPU-seconds)', file=sys.stderr, flush=True)

    #Generate the output file:
    write_pbs(args.output_tsv, windows, pbs, args.target_pop, args.reference_pop, args.outgroup_pop, args.window_size)

    t_writepbs = time.perf_counter_ns()
    tcpu_writepbs = time.process_time_ns()
    elapsed_writepbs = (t_writepbs - t_pbs) / 1_000_000_000
    elapsedcpu_writepbs = (tcpu_writepbs - tcpu_pbs) / 1_000_000_000
    elapsed_total = (t_writepbs - t_0) / 1_000_000_000
    elapsedcpu_total = (tcpu_writepbs - tcpu_0) / 1_000_000_000

    print(f'Output PBS and related values for each window to {args.output_tsv} in {elapsed_writepbs} seconds ({elapsedcpu_writepbs} CPU-seconds)', file=sys.stderr, flush=True)
    print(f'Total PBS run took {elapsed_total} seconds ({elapsedcpu_total} CPU-seconds)', file=sys.stderr, flush=True)

if __name__=="__main__":
    args = parse_args()
    run_pbs(args)

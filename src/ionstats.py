import argparse
import yaml
import os
from pathlib import Path
from preprocess import IonStatsPreprocess
from collect import IonStatsCollector
from tester import IonStatsTester
from motifs import IonStatsMotifDiscovery
from interruptions import IonStatsReadInterruptions
from summary_reader import IonStatsSummaryReader
from kmer_stats import IonStatsKmerStats
from helpers import *


def check_output_paths(out_path, step = 'collect'):
    # Create directories in out path
    # out_path/
        # /values/
            # data_type/
        # /tests/
            # data_type/
        # /significant_sequences/
            # data_type/
    make_sure_path_exists(out_path)
    make_sure_path_exists(os.path.join(out_path, "signal_alignment"))
    if step == 'collect':
        mid_dir = 'values'
    elif step == 'summary':
        mid_dir = 'read_level_values'
    elif step == 'test':
        mid_dir = 'tests'
    elif step == 'kmers':
        mid_dir = 'significant_kmers'
    elif step == 'streme':
        mid_dir = 'streme_outs'
    elif step == 'interruptions':
        mid_dir = 'interrupted_reads'
    elif step == 'stats':
        mid_dir = 'kmer_stats'
    else:
        return
    make_sure_path_exists(os.path.join(out_path, mid_dir))


def check_input_path(pth):
    base_dir = Path(__file__).resolve().parent.parent
    if pth is None:
        return None
    p = Path(pth).expanduser()
    if p.exists():
        return p.resolve()

    candidate = Path(os.path.join(base_dir, pth))
    if candidate.exists():
        print(f"[INFO] Using data paths relative to the script: {base_dir}")
        return candidate.resolve()
    
    raise FileNotFoundError(f"Could not find path: {pth}")


def run_ionstats(args):
    data_types = args.data_types
    control_id = args.control
    treated_id = args.treatment
    out_path = args.out_path
    run_id = args.run_id
    tests = args.tests
    data_path = os.path.join(out_path, "signal_alignment")

    check_output_paths(out_path)
    ref_path = check_input_path(args.ref_path)
    summary_path = check_input_path(args.summary_path)
    bam_dir = check_input_path(args.bam_path)
    read_dir = check_input_path(args.read_path)
    uc4_table_path = check_input_path(args.uc4_table_path)
    alphabet_path = check_input_path(args.alphabet)
    control_sequence_path = check_input_path(args.control_sequence_path)
    sample_ids = [control_id, treated_id]
    # The calibration control is collected like any other sample, so that the tester can
    # run the null comparison between the two controls
    if args.calibration_control:
        sample_ids.append(args.calibration_control)
    # the same identifier may be given twice (e.g. to collect one sample on its own)
    sample_ids = list(dict.fromkeys(sample_ids))
    analyses = args.analyses
    # ref_name only restricts the contigs that are used, ref_pos selects whether the
    # values are grouped by reference position instead of k-mer context
    use_refpos = args.ref_pos or (args.ref_name is not None)
    if use_refpos and 'sig_shifted_dt' in data_types:
        # The shift groups dwell times by the sequence at the motor protein, which is a
        # k-mer mode feature. Grouping by reference position already puts every dwell time
        # at the position it was measured at, so the shifted values are excluded here.
        print("Note: sig_shifted_dt is excluded in ref_pos mode, where the dwell times",
              "are already grouped by the position they were measured at. Use sig_dt.")
        data_types = [dt for dt in data_types if dt != 'sig_shifted_dt']

    if 'preprocess' in analyses:
        for si in sample_ids:
            bam_path = get_file_path(si, bam_dir)
            read_path = get_file_path(si, read_dir)
            preprocess = IonStatsPreprocess(bam_path = bam_path,
                                            read_path = read_path,
                                            ref_path = ref_path,
                                            sample_id = si,
                                            run_id = run_id,
                                            out_path = data_path,
                                            verb = args.verb)
            preprocess.signal_alignment()
    if 'collect' in analyses:
        for si in sample_ids:
            for dt in data_types:
                bam_path = get_file_path(si, data_path)
                collector = IonStatsCollector(bam_path = bam_path,
                                              sample_id = si,
                                              ref_path = ref_path,
                                              out_path = os.path.join(out_path, 'values'),
                                              run_id = run_id,
                                              data_type = dt,
                                              strand = args.strand,
                                              k = args.k,
                                              ref_name = args.ref_name,
                                              ref_pos = args.ref_pos,
                                              shift = args.shift,
                                              sensing_offset = args.sensing_offset,
                                              n_vals = args.n_vals,
                                              n_threads = args.n_threads,
                                              uc4_table_path = uc4_table_path,
                                              n_bins = args.n_bins,
                                              normalization = args.normalization,
                                              norm_table = (args.norm_table or None),
                                              norm_quantiles = args.norm_quantiles,
                                              norm_const = args.norm_const,
                                              save_every = args.save_every,
                                              verb = args.verb)
                collector.collect_values()
    if 'read_level_values' in analyses:
        check_output_paths(out_path, step = 'summary')
        for si in sample_ids:
            summary_reader = IonStatsSummaryReader(summary_path = summary_path,
                                                   sample_id = si,
                                                   out_path = os.path.join(out_path, 'read_level_values'),
                                                   run_id = run_id,
                                                   extract_columns = args.extract_columns,
                                                   id_column = args.id_column)
            summary_reader.read_summary()
    if 'test' in analyses:
        check_output_paths(out_path, step = 'test')
        check_output_paths(out_path, step = 'kmers')
        for dt in data_types:
            try:
                cv = getattr(args, f"conversion_vals.{dt}")
            except AttributeError:
                cv = None
            if dt in ['AEAD_mean', 'AEAD_std'] or use_refpos:
                it = 'dict'
            else:
                it = args.input_type
            tester = IonStatsTester(data_path = os.path.join(out_path, 'values'),
                                    test_out_path = os.path.join(out_path, 'tests'),
                                    kmer_out_path = os.path.join(out_path, 'significant_kmers'),
                                    data_type = dt,
                                    run_id = run_id,
                                    control_id = control_id,
                                    treated_id = treated_id,
                                    test_type = tests,
                                    ref_path = ref_path,
                                    ref_name = args.ref_name,
                                    ref_pos = args.ref_pos,
                                    sensing_offset = args.sensing_offset,
                                    input_type = it,
                                    strand = args.strand,
                                    k = args.k,
                                    jitter = args.jitter,
                                    js = args.js,
                                    seed = args.seed,
                                    n_splits = args.n_splits,
                                    conversion_vals = cv,
                                    alpha = args.fdr_alpha,
                                    min_run_length = args.min_run_length,
                                    calibration_id = args.calibration_control,
                                    calibration_fdr = args.calibration_fdr,
                                    peak_separation = args.peak_separation,
                                    peak_merge_distance = args.peak_merge_distance,
                                    target_bases = args.target_bases,
                                    treated_scale = args.treated_scale,
                                    treated_offset = args.treated_offset,
                                    verb = args.verb)
            tester.compare_samples()
    if 'stats' in analyses:
        check_output_paths(out_path, step = 'stats')
        for si in sample_ids:
            for dt in data_types:
                try:
                    cv = getattr(args, f"conversion_vals.{dt}")
                except AttributeError:
                    cv = None
                if dt in ['AEAD_mean', 'AEAD_std'] or use_refpos:
                    it = 'dict'
                else:
                    it = args.input_type
                kmer_stats = IonStatsKmerStats(data_path = os.path.join(out_path, 'values'),
                                               stat_out_path = os.path.join(out_path, 'kmer_stats'),
                                               data_type = dt,
                                               run_id = run_id,
                                               sample_id = si,
                                               k = args.k,
                                               input_type = it,
                                               strand = args.strand,
                                               quantiles = args.quantiles,
                                               conversion_vals = cv,
                                               n_splits = args.n_splits,
                                               ref_path = ref_path,
                                               ref_name = args.ref_name,
                                               ref_pos = args.ref_pos,
                                               sensing_offset = args.sensing_offset,
                                               verb = args.verb)
                kmer_stats.collect_values_and_stats()
    if 'motif_discovery' in analyses:
        check_output_paths(out_path, step = 'streme')
        for dt in data_types:
            for tt in tests:
                motif_discovery = IonStatsMotifDiscovery(run_id = run_id,
                                                         data_type = dt,
                                                         test_type = tt,
                                                         sequence_dir = os.path.join(out_path, 'significant_kmers'),
                                                         out_dir = os.path.join(out_path, 'streme_outs'),
                                                         alphabet = alphabet_path,
                                                         k = args.k,
                                                         control_path = control_sequence_path,
                                                         overwrite = args.overwrite,
                                                         verb = args.verb)
                motif_discovery.motif_discovery()
    if 'interruptions' in analyses:
        check_output_paths(out_path, step = 'interruptions')
        # sample_id:barcode pairs, for when the bam name and the summary barcode differ
        summary_ids = dict(p.split(":", 1) for p in args.summary_ids) if args.summary_ids else {}
        for si in sample_ids:
            bam_path = get_file_path(si, data_path)
            interruptions = IonStatsReadInterruptions(bam_path = bam_path,
                                                      summary_path = summary_path,
                                                      sample_id = si,
                                                      summary_id = summary_ids.get(si),
                                                      read_id_list = args.read_id_list,
                                                      skip_edges = args.skip_edges,
                                                      use_padding = args.use_padding,
                                                      invert_filter = args.invert_filter,
                                                      ref_path = ref_path,
                                                      out_path = os.path.join(out_path, 'interrupted_reads'),
                                                      run_id = run_id)
            interruptions.read_interruptions()
    return


def add_args_from_config(parser, config, prefix = ""):
    for key, val in config.items():
        arg_name = f"{prefix}{key}"
        if isinstance(val, dict):
            add_args_from_config(parser, val, prefix=f"{arg_name}.")
        else:
            arg_flag = f"--{arg_name}"
            arg_type = infer_type(val)
            if arg_type == 'bool':
                parser.add_argument(arg_flag, action = "store_true" if not val else "store_false",
                                    help = f"Toggle for {arg_name} (default: {val})")
            elif arg_type == 'list':
                el_type = type(val[0]) if val else str
                parser.add_argument(arg_flag, type = el_type, nargs = '+',
                                    default = val, help = f"List of {el_type.__name__} (default: {val})")
            else:
                parser.add_argument(arg_flag, type = arg_type, default = val,
                                    help = f"Default: {val}")


def parse_args(config_path = None):
    config = load_config(config_path)

    parser = argparse.ArgumentParser(description = "IonStats for analysing and comparing nanopore sequencing data." 
                                     "Check README.md for more detailed help message.")
    parser.add_argument("--config", type = str, help = "YAML config file")

    add_args_from_config(parser, config)
    
    args = parser.parse_args()
    return args


def main():
    # First, parse custom or default config, if custom not provided
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config", type=str)
    config_args, _ = config_parser.parse_known_args()

    # Next, parse args from command line
    args = parse_args(config_args.config)
    print(vars(args))
    run_ionstats(args)


if __name__ == "__main__":
    main()
    exit()

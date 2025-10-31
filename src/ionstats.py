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
    if step == 'summary':
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
    # Which of the input should be as arguments, and which as config
    # Maybe the default use should be that a path should be provided for 
    # bam and pod5 and one identifier for control and one for treated.
    # For more sophisticated uses (multiple inputs per treatment or multiple treatments)
    # the identifiers could be set up in the config.
    # The config should be something that a normal user should not have to touch or even specify.
    analyses = args.analyses

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
                                              n_vals = args.n_vals,
                                              n_threads = args.n_threads,
                                              uc4_table_path = uc4_table_path,
                                              n_bins = args.n_bins,
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
            tester = IonStatsTester(data_path = os.path.join(out_path, 'values'),
                                    test_out_path = os.path.join(out_path, 'tests'),
                                    kmer_out_path = os.path.join(out_path, 'significant_kmers'),
                                    data_type = dt,
                                    run_id = run_id,
                                    control_id = control_id,
                                    treated_id = treated_id,
                                    test_type = tests,
                                    input_type = args.input_type,
                                    strand = args.strand,
                                    k = args.k,
                                    jitter = args.jitter,
                                    js = args.js,
                                    n_splits = args.n_splits,
                                    conversion_vals = args.conversion_vals,
                                    verb = args.verb)
            tester.compare_samples()
    if 'stats' in analyses:
        check_output_paths(out_path, step = 'stats')
        for si in sample_ids:
            for dt in data_types:
                kmer_stats = IonStatsKmerStats(data_path = os.path.join(out_path, 'values'),
                                               stat_out_path = os.path.join(out_path, 'kmer_stats'),
                                               data_type = dt,
                                               run_id = run_id,
                                               sample_id = si,
                                               k = args.k,
                                               input_type = args.input_type,
                                               strand = args.strand,
                                               quantiles = args.quantiles,
                                               conversion_vals = args.conversion_vals,
                                               n_splits = args.n_splits,
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
        for si in sample_ids:
            bam_path = get_file_path(si, data_path)
            interruptions = IonStatsReadInterruptions(bam_path = bam_path,
                                                      summary_path = summary_path,
                                                      sample_id = si,
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

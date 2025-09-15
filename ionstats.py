import argparse
import yaml
import os
from preprocess import IonStatsPreprocess
from collect import IonStatsCollector
from test import IonStatsTester
from motif_discovery import IonStatsMotifDiscovery
from helpers import *


def check_output_paths(out_path, data_types, step = 'collect'):
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
    elif step == 'test':
        mid_dir = 'tests'
    elif step == 'kmers':
        mid_dir = 'significant_kmers'
    elif step == 'streme':
        mid_dir = 'streme_outs'
    make_sure_path_exists(os.path.join(out_path, mid_dir))
    for dt in data_types:
        new_path = os.path.join(out_path, mid_dir, dt)
        make_sure_path_exists(new_path)
    pass


def run_ionstats(args):
    config_file = args.config
    with open(config_file) as file:
        config = yaml.safe_load(file)

    data_types = args.data_types.split(",")
    control_id = args.control
    treated_id = args.treatment
    out_path = args.out_path
    ref_path = args.ref_path
    run_id = args.run_id
    tests = args.tests
    bam_dir = args.bam_path
    read_dir = args.read_path
    data_path = os.path.join(out_path, "signal_alignment")
    check_output_paths(out_path, data_types)

    sample_ids = [control_id, treated_id]
    # Which of the input should be as arguments, and which as config
    # Maybe the default use should be that a path should be provided for 
    # bam and pod5 and one identifier for control and one for treated.
    # For more sophisticated uses (multiple inputs per treatment or multiple treatments)
    # the identifiers could be set up in the config.
    # The config should be something that a normal user should not have to touch or even specify.
    analyses = config["analyses"]

    if 'preprocess' in analyses:
        for si in sample_ids:
            bam_path = get_file_path(si, bam_dir)
            read_path = get_file_path(si, read_dir)
            print(bam_path, read_path)
            print(bam_dir, read_dir)
            preprocess = IonStatsPreprocess(bam_path = bam_path,
                                            read_path = read_path,
                                            ref_path = ref_path,
                                            sample_id = si,
                                            run_id = run_id,
                                            out_path = data_path,
                                            verb = config["verb"])
            preprocess.signal_alignment()
    if 'collect' in analyses:
        for si in sample_ids:
            for dt in data_types:
                bam_path = get_file_path(si, data_path)
                collector = IonStatsCollector(bam_path = bam_path,
                                              sample_id = si,
                                              ref_path = ref_path,
                                              out_path = os.path.join(out_path, 'values', dt),
                                              run_id = run_id,
                                              data_type = dt,
                                              strand = config["strand"],
                                              k = config["k"],
                                              ref_name = config["ref_name"],
                                              n_vals = config["n_vals"],
                                              n_threads = config["n_threads"],
                                              uc4_table_path = config["uc4_table_path"],
                                              n_bins = config["n_bins"],
                                              save_every = config["save_every"],
                                              verb = config["verb"])
                collector.collect_values()
    check_output_paths(out_path, data_types, step = 'test')
    check_output_paths(out_path, data_types, step = 'kmers')
    if 'test' in analyses:
        for dt in data_types:
            tester = IonStatsTester(data_path = os.path.join(out_path, 'values', dt),
                                    test_out_path = os.path.join(out_path, 'tests', dt),
                                    kmer_out_path = os.path.join(out_path, 'significant_kmers', dt),
                                    data_type = dt,
                                    run_id = run_id,
                                    control_id = control_id,
                                    treated_id = treated_id,
                                    test_type = tests,
                                    input_type = config["input_type"],
                                    strand = config["strand"],
                                    k = config["k"],
                                    jitter = config["jitter"],
                                    js = config["js"],
                                    n_splits = config["n_splits"],
                                    conversion_vals = config["conversion_vals"],
                                    verb = config["verb"])
            tester.compare_samples()
    check_output_paths(out_path, data_types, step = 'streme')
    if 'motif_discovery' in analyses:
        for tt in tests.split(","):
            for dt in data_types:
                motif_discovery = IonStatsMotifDiscovery(run_id=run_id,
                                                        sequence_dir=os.path.join(out_path, 'significant_kmers', dt),
                                                        out_dir=os.path.join(out_path, 'streme_outs', dt),
                                                        alphabet=config["alphabet"],
                                                        test_type = tt,
                                                        k =config["k"],
                                                        overwrite=config["overwrite"],
                                                        verb = config["verb"])
            motif_discovery.motif_discovery()
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--control", action="store", type=str, required=True)
    parser.add_argument("--treatment", action="store", type=str, required=True)
    parser.add_argument("--bam_path", action="store", type=str, required=True)
    parser.add_argument("--read_path", action="store", type=str, required=True)
    parser.add_argument("--out_path", action="store", type=str, required=True)
    parser.add_argument("--ref_path", action="store", type=str, required=True)
    parser.add_argument("--run_id", action="store", type=str, required=True)
    parser.add_argument("--tests", action="store", type=str, default="ks,cvm")
    parser.add_argument("-c", "--config", action="store", type=str, default="default_config.yaml")
    parser.add_argument("--data_types", action="store", type=str, default="qs,sig_mean,sig_std,sig_dt")
    args = parser.parse_args()
    run_ionstats(args)


if __name__ == "__main__":
    main()
    exit()

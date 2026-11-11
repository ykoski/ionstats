import numpy as np
import os
import pickle
import pandas as pd

from helpers import *

class IonStatsKmerStats():

    def __init__(self,
                 data_path,
                 stat_out_path,
                 data_type,
                 run_id,
                 sample_id,
                 k,
                 input_type = 'dict',
                 strand = 'both',
                 quantiles = True,
                 conversion_vals = [-5., 0.00015259021896696422],
                 n_splits = 1,
                 ref_path = None,
                 ref_name = None,
                 ref_pos = False,
                 sensing_offset = 2,
                 verb = True):
        self.data_path = data_path
        self.data_type = data_type
        self.out_path = stat_out_path
        self.run_id = run_id
        self.sample_ids = sample_id.split(",") # Implement functionality for sample groups later
        self.k = k
        self.input_type = input_type
        self.strand = strand
        self.quantiles = quantiles
        self.n_splits = n_splits
        self.ref_path = ref_path
        self.ref_name = ref_name
        # ref_name only restricts which contigs are used, ref_pos selects the mode
        self.ref_names = None
        if self.ref_name is not None:
            self.ref_names = set(self.ref_name.split(","))
        self.use_refpos = ref_pos or (self.ref_names is not None)
        self.sensing_offset = sensing_offset
        if self.use_refpos:
            # Splitting is done over the k-mer set, which is not used in ref_pos mode
            self.n_splits = 1
        self.conversion_vals = conversion_vals
        self.verb = verb
        if self.verb:
            print(f"Set up k-mer stats for samples {self.sample_ids} and for data type {self.data_type}.")


    def load_distribution(self, ids, kmer_set):
        # membership is tested once per k-mer per file, so make it O(1)
        if kmer_set is not None and not isinstance(kmer_set, (set, frozenset)):
            kmer_set = set(kmer_set)
        files = get_matching_files(self.data_path, self.strand, self.data_type, ids)
        print(files)
        data = None
        for f in files:
            if f.endswith('.npy'):
                if self.input_type == 'dict':
                    print("Warning! The input type does not match loaded data!")
                new_data = np.load(f)
                if data is None:
                    data = new_data
                else:
                    data += new_data
            elif f.endswith('.pickle'):
                if self.input_type == 'array':
                    print("Warning! The input type does not match loaded data!")
                with open(f, 'rb') as handle:
                    new_data = pickle.load(handle)
                if data is None:
                    data = new_data if kmer_set is None else \
                           {k: v for k, v in new_data.items() if k in kmer_set}
                else:
                    data = update_dict(data, new_data, kmer_set = kmer_set, np_arrays=True)
        return data
    

    def compute_stat_df(self, data):
        means, stds, medians, index, nvals = [], [], [], [], []
        lq1s, hq1s, lq5s, hq5s = [], [], [], []
        contigs, positions, ref_kmer_list = [], [], []
        if self.input_type == 'dict':
            # In ref_pos mode the keys are (contig, ref_pos) tuples, which sort by
            # contig and then by position
            iter_set = sorted(data.keys())
            if self.use_refpos and self.ref_names is not None:
                iter_set = [key for key in iter_set if key[0] in self.ref_names]
            index_name = 'kmer'
        elif self.input_type == 'array':
            iter_set = range(len(data))
            index_name = 'ref_pos'
        if self.use_refpos and self.ref_path is not None:
            ref_kmers = RefKmerFetcher(self.ref_path, self.k, self.sensing_offset)
        else:
            ref_kmers = None
        for i in iter_set:
            if self.input_type == 'dict':
                values = np.array(data[i], dtype = np.float64)
            elif self.input_type == 'array':
                a = np.array(data[i])
                values = counts_to_array(a, self.conversion_vals)
            if len(values) <= 1:
                continue
            means.append(np.mean(values))
            stds.append(np.std(values))
            medians.append(np.median(values))
            if self.quantiles:
                if len(values) < 100:
                    print(f"Warning! Less than 100 observations for kmer: {i} in {self.sample_ids}."
                          "You should be careful when interpreting the quantile values.")
                q1 = np.quantile(values, (0.01,0.99))
                q5 = np.quantile(values, (0.05,0.95))
                lq1s.append(q1[0])
                lq5s.append(q5[0])
                hq1s.append(q1[1])
                hq5s.append(q5[1])
            nvals.append(len(values))
            if self.use_refpos and self.input_type == 'dict':
                contig, pos = i
                contigs.append(contig)
                positions.append(pos)
                if ref_kmers is not None:
                    ref_kmer_list.append(ref_kmers.kmer_at(contig, pos))
            else:
                index.append(i)
        if self.use_refpos and self.input_type == 'dict':
            df_data = {'contig' : contigs, 'ref_pos' : positions}
            if ref_kmers is not None:
                df_data['kmer'] = ref_kmer_list
            sort_cols = ['contig', 'ref_pos']
        else:
            df_data = {index_name : index}
            sort_cols = [index_name]
        df_data.update({'mean' : means, 'std' : stds, 'median' : medians})
        if self.quantiles:
            df_data.update({'lq1' : lq1s, 'hq1' : hq1s, 'lq5' : lq5s, 'hq5' : hq5s})
        df_data['nvals'] = nvals
        df = pd.DataFrame(data = df_data)
        df = df.sort_values(by=sort_cols)
        return df


    def save_stat_dataframe(self, df, sid = None):
        if sid is None:
            path_list = [self.sample_ids, self.run_id, self.data_type]
        else:
            path_list = [sid, self.run_id, self.data_type]
        if self.strand != 'both':
            path_list.append(self.strand)
        fname = os.path.join(self.out_path, "_".join(path_list)) + ".csv"
        df.to_csv(fname)
        if self.verb:
            print(f"Saved the stat output to {fname}")


    def collect_values_and_stats(self):
        if self.n_splits > 1:
            kmer_set = get_all_kmers(k=self.k)
            split_size  = len(kmer_set) // self.n_splits

        for sid in self.sample_ids:
            for split in range(self.n_splits):
                if self.n_splits > 1:
                    kmers = kmer_set[split_size*split:split_size*(split+1)]
                else:
                    kmers = None
                data = self.load_distribution([sid], kmers)
                # Compute stats from data
                new_df = self.compute_stat_df(data)
                if split >= 1:
                    df = pd.concat([df, new_df])
                else:
                    df = new_df
            self.save_stat_dataframe(df, sid)
            if self.verb:
                print(f"Finished stat computation of {sid}")

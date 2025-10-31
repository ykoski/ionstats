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
        self.conversion_vals = conversion_vals
        self.verb = verb
        if self.verb:
            print(f"Set up k-mer stats for samples {self.sample_id} and for data type {self.data_type}.")


    def load_distribution(self, ids, kmer_set):
        #files = self.get_matching_files(ids)
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
                    data = new_data
                else:
                    data = update_dict(data, new_data, kmer_set = kmer_set, np_arrays=True)
        return data
    

    def compute_stat_df(self, data):
        means, stds, medians, index, nvals = [], [], [], [], []
        if self.quantiles:
            lq1s, hq1s, lq5s, hq5s = [], [], [], []
        if self.input_type == 'dict':
            iter_set = data.keys()
            index_name = 'kmer'
        elif self.input_type == 'array':
            iter_set = range(len(data))
            index_name = 'ref_pos'
        for i in iter_set:
            if self.input_type == 'dict':
                values = np.array(data[i], dtype = np.float64)
            elif self.input_type == 'array':
                a = np.array(data[i])
                values = counts_to_array(a, self.conversion_vals)
            means.append(np.mean(values))
            stds.append(np.std(values))
            medians.append(np.median(values))
            if self.quantiles:
                if len(values) < 100:
                    print(f"Warning! Less than 100 observations for kmer: {i} in {self.sample_id}." 
                          "You should be careful when interpreting the quantile values.")
                q1 = np.quantile(values, (0.01,0.99))
                q5 = np.quantile(values, (0.05,0.95))
                lq1s.append(q1[0])
                lq5s.append(q5[0])
                hq1s.append(q1[1])
                hq5s.append(q5[1])
            nvals.append(len(values))
            index.append(i)
        df = pd.DataFrame(data = {index_name : index,
                                'mean' : means,
                                'std' : stds,
                                'median' : medians,
                                'lq1' : lq1s,
                                'hq1' : hq1s,
                                'lq5' : lq5s,
                                'hq5' : hq5s,
                                'nvals' : nvals})
        if self.quantiles:
            df = df.assign(lq1 = lq1s, hq1 = hq1s, lq5 = lq5s, hq5 = hq5s)
        df = df.sort_values(by=index_name)
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
        kmer_set = get_all_kmers(k=self.k)
        if self.n_splits > 1:
            split_size  = len(kmer_set) // self.n_splits

        for sid in self.sample_ids:
            for split in range(self.n_splits):
                if self.n_splits > 1:
                    kmers = kmer_set[split_size*split:split_size*(split+1)]
                else:
                    kmers = None
                data = self.load_distribution(sid, kmers)
                # Compute stats from data
                new_df = self.compute_stat_df(data)
                if split >= 1:
                    df = pd.concat(df, new_df)
                else:
                    df = new_df
            self.save_stat_dataframe(df, sid)
            if self.verb:
                print(f"Finished stat computation of {sid}")

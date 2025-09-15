import numpy as np
import scipy.stats as scp
import pickle
import os
import pandas as pd
import pysam

from statsmodels.stats.multitest import fdrcorrection

from helpers import *

class IonStatsTester():

    def __init__(self,
                 data_path,
                 test_out_path,
                 kmer_out_path,
                 data_type,
                 run_id,
                 control_id,
                 treated_id,
                 test_type,
                 input_type = 'dict',
                 strand = 'both',
                 k = 9,
                 jitter = False,
                 js = 0.01,
                 n_splits = 1,
                 conversion_vals = [-5., 0.00015259021896696422],
                 alpha = 0.05,
                 verb = True):
        self.data_path = data_path
        self.test_out_path = test_out_path
        self.kmer_out_path = kmer_out_path
        self.data_type = data_type
        self.run_id = run_id
        self.control_ids = control_id.split(",")
        self.treated_ids = treated_id.split(",")
        self.input_type = input_type
        self.test_types = test_type.split(",")
        self.strand = strand
        self.jitter = jitter
        self.js = js
        self.alpha = alpha
        self.n_splits = n_splits
        if self.input_type == 'array':
            self.n_splits = 1
        self.k = k
        self.conversion_vals = conversion_vals
        self.verb = verb
        if self.verb:
            print(f"Set up tester for samples {control_id, treated_id} and for data type {self.data_type} using {test_type} test.")


    def get_matching_files(self, ids):
        matching_files = []
        all_files = os.listdir(self.data_path)
        for id in ids:
            for f in all_files:
                full_path = os.path.join(self.data_path, f)
                if self.strand != 'both':
                    if id in f and self.strand in f and self.data_type in f:
                        matching_files.append(full_path)
                else:
                    if id in f and self.data_type in f:
                        matching_files.append(full_path)
        return matching_files


    def load_distribution(self, ids, kmer_set):
        files = self.get_matching_files(ids)
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


    def test_samples(self, control_data, treated_data, test_type):
        if self.input_type == 'dict':
            keys = list(control_data.keys() & treated_data.keys())
            n_values = len(keys)
            id_col = 'kmer'
        elif self.input_type == 'array':
            n_values = control_data.shape[0]
            id_col = 'refpos'
        test_values, p_values, ids = [], [], []
        for i in range(n_values):
            if self.input_type == 'dict':
                key = keys[i]
                c_arr = control_data[key]
                t_arr = treated_data[key]
                id = key
            elif self.input_type == 'array':
                c_arr = counts_to_array(control_data[i], self.conversion_vals)
                t_arr = counts_to_array(treated_data[i], self.conversion_vals)
                id = i
            if self.jitter:
                t_arr = t_arr.astype(float) + np.random.uniform(-self.js, self.js, size=t_arr.shape)
                c_arr = c_arr.astype(float) + np.random.uniform(-self.js, self.js, size=c_arr.shape)
            if len(t_arr) <= 1 or len(c_arr) <= 1:
                continue
            if test_type == 'ks':
                stat = scp.ks_2samp(t_arr, c_arr)
            if test_type == 'cvm':
                stat = scp.cramervonmises_2samp(t_arr, c_arr)
            test_values.append(stat.statistic)
            p_values.append(stat.pvalue)
            ids.append(id)
        df = pd.DataFrame(data = {id_col : ids,
                                  test_type : test_values,
                                  'pvalue' : p_values})
        df = df.sort_values(by = id_col)
        return df


    def save_test_dataframe(self, df, test):
        path_list = [self.run_id, self.data_type, test]
        if self.strand != 'both':
            path_list.append(self.strand)
        fname = os.path.join(self.test_out_path, "_".join(path_list)) + ".csv"
        df.to_csv(fname)
        if self.verb:
            print(f"Saved the test output to {fname}")


    def write_significant_kmers(self, df, test):
        filt, _ = fdrcorrection(df['pvalue'], alpha = self.alpha)
        df_f = df[filt]
        write_out_sequences(list(df_f.kmer), 
                            os.path.join(self.kmer_out_path, "_".join([self.run_id, self.data_type, test]) + '.fasta'))


    def compare_samples(self):
        kmer_set = get_all_kmers(k=self.k)
        if self.n_splits > 1:
            split_size  = len(kmer_set) // self.n_splits
        results = {}
        for split in range(self.n_splits):
            if self.n_splits > 1:
                kmers = kmer_set[split_size*split:split_size*(split+1)]
            else:
                kmers = None
            control_data = self.load_distribution(self.control_ids, kmers)
            treated_data = self.load_distribution(self.treated_ids, kmers)
            for tt in self.test_types:
                test_df = self.test_samples(control_data, treated_data, tt)
                if split >= 1:
                    results[tt] = pd.concat(results[tt], test_df)
                else:
                    results[tt] = test_df
        for test, df in results.items():
            self.save_test_dataframe(df, test)
            if self.input_type == 'dict':
                self.write_significant_kmers(df, test)
        if self.verb:
            print(f"Finished testing of {self.control_ids} and {self.treated_ids} using {self.test_types}")


                

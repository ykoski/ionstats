import numpy as np
import pysam
import os
import pickle

from helpers import *

class IonStatsCollector():

    def __init__(self,
                 bam_path,
                 sample_id,
                 ref_path,
                 out_path,
                 run_id,
                 data_type,
                 strand = 'both',
                 k = 9,
                 ref_name = None,
                 n_vals = 65536,
                 n_threads = 1,
                 uc4_table_path = None,
                 n_bins = None,
                 save_every = None,
                 verb = True):
        self.bam_path = bam_path
        self.sample_id = sample_id
        self.ref_path = ref_path
        self.alignment = pysam.AlignmentFile(self.bam_path, threads = n_threads)
        self.reference = pysam.FastaFile(self.ref_path)
        self.reflens = self.reference.lengths
        self.len_dict = dict(zip(self.reference.references, self.reference.lengths))
        self.out_path = out_path
        self.run_id = run_id
        self.data_type = data_type
        self.ref_name = ref_name
        self.uc4_table_path = uc4_table_path
        self.strand = strand
        self.use_refpos = False
        self.k = k
        self.h = self.k // 2
        if self.ref_name is not None:
            self.use_refpos = True
        self.n_vals = n_vals
        self.signal_coeff = 10 / (self.n_vals-1)
        self.save_every = save_every
        self.n_bins = n_bins
        self.binning = False
        if self.n_bins is not None:
            self.binning = True
            kmers = get_all_kmers(self.k)
            self.kmer_map = {k : i for i, k in enumerate(kmers)}
            self.bw = n_vals / self.n_bins # This can be done since the range of values is 
            # from 0 to n_vals
        if self.uc4_table_path is not None:
            if self.data_type != 'AEAD_mean' and self.data_type != 'AEAD_std':
                print(f"Warning: uc4_table path passed with data_type {self.data_type}.",
                      "The kmer table is not used or loaded.")
            else:
                # Make a modification for the ref_pos case, where only the kmer values 
                # of different reference positions are fetched
                mean_dict, std_dict = get_kmer_mean_and_std(self.uc4_table_path)
                self.uc4_means = mean_dict
                self.uc4_stds = std_dict
        self.verb = verb
        if self.verb:
            print(f"Set up collector for sample {self.sample_id} and for data type {self.data_type}.")


    def init_datastructure(self):
        if self.binning:
            n_vals = self.n_bins
            if not self.use_refpos:
                return np.zeros((self.k**4, n_vals))
        else:
            n_vals = self.n_vals
        if not self.use_refpos:
            return {}
        else:
            reflen = self.len_dict[self.refname]
            return np.zeros((reflen, n_vals))


    def get_read_values(self, read, read_mapping):
        if self.data_type == 'qs':
            refname = self.alignment.get_reference_name(read.reference_id)
            pairs = read.get_aligned_pairs(matches_only = True)
            al_start = read.reference_start
            al_end = read.reference_end
            qs = read.get_forward_qualities()
            qs_seq_start = al_start - self.h
            qs_seq_end = al_end + self.h
            start_pad, end_pad = 0, 0
            if qs_seq_start < 0:
                start_pad = abs(qs_seq_start)
                qs_seq_start = 0
            contig_len = self.len_dict[refname]
            if qs_seq_end > contig_len:
                end_pad = abs(contig_len - qs_seq_end)
                qs_seq_end = contig_len
            # Fetch new sequence that is always in forward orientation and 
            # is adjusted for qs_seq_start and qs_seq_end
            qs_sequence = self.reference.fetch(refname, qs_seq_start, qs_seq_end)
            qs_kmers = get_kmers(qs_sequence, 9)
            qs = qs[start_pad:len(qs)-end_pad]
            values, read_mapping = self.get_matched_quality_scores(pairs, qs, qs_kmers, qs_seq_start)
        elif self.data_type == "sig_mean" or self.data_type == 'AEAD_mean':
            values = np.array(read.get_tag('uc'), dtype = int)
        elif self.data_type == 'sig_std' or self.data_type == 'AEAD_std':
            values = np.array(read.get_tag('ud'), dtype = int)
        elif self.data_type == 'sig_dt' or self.data_type == 'sig_shifted_dt':
            ul = np.array(read.get_tag('ul'), dtype=int)
            dts = reformat_ul(ul)
            # Get rid of trimmed signal samples
            values = dts[1:]
            values = values[values >= 0]
        else:
            print("Unknown data_type passed. Exiting...")
            exit()
        values = self.kmer_idx_to_amps(values)
        if self.data_type == 'AEAD_mean':
            exp_means, _ = get_expected_values_for_read(read_mapping, self.uc4_means)
            mean_diffs = np.abs(exp_means - values)
            values = compute_AEAD_values(mean_diffs, k = self.k)
        elif self.data_type == 'AEAD_std':
            _, exp_stds = get_expected_values_for_read(read_mapping, self.uc4_stds)
            std_diffs = np.abs(exp_stds - values)
            values = compute_AEAD_values(std_diffs, k = self.k)
        return values, read_mapping


    def get_matched_quality_scores(self, pairs, quals, kmers, ref_window_start, max_qs=None):
        values = []
        out_kmers = []
        for (sp, rp), q in zip(pairs, quals):
            kmer_idx = rp - ref_window_start
            if 0 <= kmer_idx < len(kmers):
                if max_qs is not None and q > max_qs:
                    q = max_qs
                values.append(q)
                out_kmers.append(kmers[kmer_idx])
        return values, out_kmers


    def get_kmers_or_refpos(self, read):
        ur = read.get_tag('ur')
        rev = read.is_reverse
        refname = self.alignment.get_reference_name(read.reference_id)
        if not self.use_refpos:
            sequence = self.reference.fetch(refname, ur[0], ur[1])
            if rev:
                sequence = reverse_complement(sequence)
            kmers = get_kmers(sequence, self.k)
            return kmers
        else:
            refpos = np.arange(ur[0]+self.h, ur[1]-self.h)
            if rev:
                refpos = refpos[::-1]
            return refpos
        
    
    def kmer_idx_to_amps(self, values):
        if self.binning:
            values = values + self.n_vals // 2
            values = values // self.bw
        if self.use_refpos:
            if self.data_type != 'AEAD_mean' and self.data_type != 'AEAD_std':
                return values    
        elif self.data_type == 'sig_dt' or self.data_type == 'sig_shifted_dt' or self.data_type == 'qs':
            return values
        values = np.array(values * self.signal_coeff, dtype = np.float16)
        return values


    def collect_values(self):
        # Init data dictionary
        data = self.init_datastructure()
        n_reads = 0
        if self.save_every is not None:
            file_idx = 0
        else:
            file_idx = None
        for read in self.alignment.fetch(until_eof=True):
            # Always skip non-primary alignments
            if read.is_secondary or read.is_supplementary:
                continue
            ur = read.get_tag('ur')
            if len(ur) > 2:
                # Skip split reads
                continue
            rev = read.is_reverse
            if (self.strand == 'forward' and rev) or (self.strand == 'reverse' and not rev):
                continue
            read_mapping = self.get_kmers_or_refpos(read)
            values, read_mapping = self.get_read_values(read, read_mapping)
            # NOTE: for AEAD values, the array will be length len(read_mapping) - 8
            if self.data_type == 'AEAD_mean' or self.data_type == 'AEAD_std':
                read_mapping = read_mapping[self.h:len(read_mapping)-self.h+1]
            assert len(read_mapping) == len(values), f"Mapping should be the same length as values {len(read_mapping)} != {len(values)}"
            if not self.use_refpos:
                if self.binning:
                    for i in range(len(values)):
                        kmer_idx = self.kmer_map(read_mapping[i])
                        data[kmer_idx, values[i]] += 1
                else:
                    for i in range(len(values)):
                        add_to_dict(data, read_mapping[i], values[i])
            else:
                for i in range(len(values)):
                    data[read_mapping[i], values[i]] += 1
            n_reads += 1
            if self.save_every is not None:
                if n_reads % self.save_every == 0:
                    self.save_values(data, file_idx)
                    data = self.init_datastructure()
                    file_idx += 1
        self.save_values(data, file_idx)
        if self.verb:
            print(f"Collection of values for {self.sample_id} finished.")


    def save_values(self, data, file_idx = None):
        path_list =  [self.sample_id, self.run_id, self.data_type]
        if self.strand != "both":
            path_list.append(self.strand)
        if file_idx is not None:
            path_list.append(str(file_idx))
        fname = os.path.join(self.out_path, "_".join(path_list))
        if self.use_refpos or self.binning:
            fname += ".npy"
            np.save(fname, data)
        else:
            trim_dict(data)
            fname += ".pickle"
            with open(fname, 'wb') as handle:
                pickle.dump(data, handle, protocol = pickle.HIGHEST_PROTOCOL)
        if self.verb:
            print(f"Saved collector output to {fname}")
import numpy as np
import scipy.stats as scp
import pickle
import os
import warnings
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
                 ref_path,
                 ref_name,
                 ref_pos = False,
                 sensing_offset = 2,
                 input_type = 'dict',
                 strand = 'both',
                 k = 9,
                 jitter = False,
                 js = 0.01,
                 seed = 42,
                 n_splits = 1,
                 conversion_vals = [-5., 0.00015259021896696422],
                 alpha = 0.05,
                 min_run_length = 1,
                 calibration_id = None,
                 calibration_fdr = 0.01,
                 peak_separation = 2,
                 peak_merge_distance = 5,
                 target_bases = None,
                 treated_scale = 1.0,
                 treated_offset = 0.0,
                 verb = True):
        self.data_path = data_path
        self.test_out_path = test_out_path
        self.kmer_out_path = kmer_out_path
        self.data_type = data_type
        self.run_id = run_id
        self.control_ids = control_id.split(",")
        self.treated_ids = treated_id.split(",")
        self.input_type = input_type
        self.test_types = test_type
        self.strand = strand
        self.jitter = jitter
        self.js = js
        # Jitter breaks ties in the discrete data types so that the continuity
        # assumption of KS and CvM approximately holds. It is RANDOM, so the
        # generator is seeded and local: an unseeded global np.random made two
        # runs of the same comparison disagree, which is a worse problem than
        # the ties. Note AD with midrank is tie corrected and needs no jitter.
        self.seed = int(seed)
        self.rng = np.random.default_rng(self.seed)
        self.alpha = alpha
        self.n_splits = n_splits
        self.ref_path = ref_path
        self.ref_name = ref_name
        # ref_name only restricts which contigs are tested, ref_pos selects the mode
        self.ref_names = None
        if self.ref_name is not None:
            self.ref_names = set(self.ref_name.split(","))
        self.use_refpos = ref_pos or (self.ref_names is not None)
        # In ref_pos mode the positions name the base that dominates the signal, which
        # sits sensing_offset bases from the k-mer centre
        self.sensing_offset = sensing_offset
        if self.input_type == 'array' or self.use_refpos:
            # Splitting is done over the k-mer set, which is not used in ref_pos mode
            self.n_splits = 1
        self.k = k
        self.h = k // 2
        self.conversion_vals = conversion_vals
        # A modification perturbs every k-mer that contains it, so a real hit shows up as a
        # run of neighbouring significant positions, while library differences tend to give
        # isolated ones. Only used in ref_pos mode, where the positions are adjacent.
        self.min_run_length = int(min_run_length)
        self.run_gap = 1
        # Primary calibration: a second control sample. The same test is run between the
        # two controls, where no modification exists, and the resulting distribution of the
        # test statistic is the empirical null that sets the threshold. This works in both
        # modes, since it only needs the statistic of the same key in a null comparison.
        self.calibration_ids = calibration_id.split(",") if calibration_id else None
        # The null is counted rather than assumed: the threshold is the most permissive one
        # at which (null values above it) / (observed values above it) stays within
        # calibration_fdr. This is a lower bound on the true rate, because a comparison of
        # two controls carries no modification-induced degradation of its own.
        self.calibration_fdr = float(calibration_fdr)
        # A modification perturbs a whole window of positions, so the called positions form
        # smears rather than points. Peaks are the local maxima of the test statistic among
        # them, and a peak is primary when no stronger peak sits close by. The default
        # separation of 2 keeps neighbouring peaks one base apart, so that two modifications
        # in e.g. a CGCG context can still be reported separately.
        self.peak_separation = int(peak_separation)
        self.peak_merge_distance = int(peak_merge_distance)
        # When the modified base is known, a peak is moved onto the nearest position that
        # could carry the modification instead of being suppressed by distance.
        self.treated_scale = float(treated_scale)
        self.treated_offset = float(treated_offset)
        self.correcting = (self.treated_scale != 1.0) or (self.treated_offset != 0.0)
        self.target_bases = set(target_bases.upper()) if target_bases else None
        self.snap_distance = max(self.peak_separation - 1, 1)
        self.verb = verb
        if self.verb:
            print(f"Set up tester for samples {control_id, treated_id} and for data type {self.data_type} using {test_type} test.")
            if self.correcting:
                how = f"(x - {self.treated_offset:.6g}) / {self.treated_scale:.6g}"
                print(f"Treated values mapped onto the control scale before testing using"
                      f" {how}. The calibration pair is left uncorrected.")


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
        # membership is tested once per k-mer per file, so make it O(1)
        if kmer_set is not None and not isinstance(kmer_set, (set, frozenset)):
            kmer_set = set(kmer_set)
        #files = self.get_matching_files(ids)
        files = get_matching_files(self.data_path, self.strand, self.data_type, ids)
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
                    # The k-mer subset has to be applied here too, not only when merging.
                    # With a single chunk file this branch was the only one taken, so the
                    # full dictionary came back whatever subset was asked for, and every
                    # split then produced stats for all k-mers instead of its own slice.
                    data = new_data if kmer_set is None else \
                           {k: v for k, v in new_data.items() if k in kmer_set}
                else:
                    data = update_dict(data, new_data, kmer_set = kmer_set, np_arrays=True)
        return data


    def test_samples(self, control_data, treated_data, test_type, correct = True):
        keys = None
        if self.input_type == 'dict':
            # In ref_pos mode the keys are (contig, ref_pos) tuples, which sort by
            # contig and then by position
            keys = sorted(control_data.keys() & treated_data.keys())
            if self.use_refpos and self.ref_names is not None:
                keys = [key for key in keys if key[0] in self.ref_names]
            n_values = len(keys)
            id_col = 'kmer'
        elif self.input_type == 'array':
            n_values = control_data.shape[0]
            id_col = 'ref_pos'
        if self.use_refpos and self.input_type == 'dict':
            ref_kmers = RefKmerFetcher(self.ref_path, self.k, self.sensing_offset)
        else:
            ref_kmers = None
        test_values, p_values, ids = [], [], []
        mean_diffs, median_diffs, n_treated, n_control, scaled = [], [], [], [], []
        contigs, positions, ref_kmer_list = [], [], []
        for i in range(n_values):
            if keys is not None:
                key = keys[i]
                c_arr = control_data[key]
                t_arr = treated_data[key]
            else:
                key = i
                c_arr = counts_to_array(control_data[i], self.conversion_vals)
                t_arr = counts_to_array(treated_data[i], self.conversion_vals)
            if correct and self.correcting:
                # onto the control scale; applied before any jitter so the tie breaking
                # noise keeps the size it was configured with
                t_arr = np.asarray(t_arr, dtype = float)
                t_arr = (t_arr - self.treated_offset) / self.treated_scale
            if self.jitter:
                t_arr = t_arr.astype(float) + self.rng.uniform(-self.js, self.js, size=t_arr.shape)
                c_arr = c_arr.astype(float) + self.rng.uniform(-self.js, self.js, size=c_arr.shape)
            if len(t_arr) <= 1 or len(c_arr) <= 1:
                continue
            if test_type == 'ks':
                stat = scp.ks_2samp(t_arr, c_arr)
            elif test_type == 'cvm':
                stat = scp.cramervonmises_2samp(t_arr, c_arr)
            elif test_type == 'ad':
                # midrank = True is the tie corrected Scholz-Stephens variant, which is
                # the one that applies to discrete data: qs and the dwell times take only
                # a handful of distinct values, and the uncorrected version assumes a
                # continuous distribution. The statistic scipy returns is already
                # standardised by its own null mean and standard deviation for the actual
                # sample sizes and tie pattern, so unlike the KS statistic it needs no
                # further scaling to be comparable between k-mers of different coverage.
                try:
                    with warnings.catch_warnings():
                        # scipy warns whenever the p value falls outside its table, which
                        # is almost always here; the statistic is unaffected
                        warnings.simplefilter('ignore')
                        stat = scp.anderson_ksamp([t_arr, c_arr], midrank = True)
                except ValueError:
                    # every observation identical, so there is nothing to compare
                    continue
            else:
                raise ValueError(f"Unknown test type '{test_type}'. Use ks, cvm or ad.")
            test_values.append(stat.statistic)
            p_values.append(stat.pvalue)
            n_treated.append(len(t_arr))
            n_control.append(len(c_arr))
            if test_type == 'ks':
                # the argument of the Kolmogorov limit law: no underflow, coverage aware
                scaled.append(stat.statistic * np.sqrt(len(t_arr) * len(c_arr)
                                                       / (len(t_arr) + len(c_arr))))
            # The test statistic has no direction, so report the effect size as well
            mean_diffs.append(np.mean(t_arr) - np.mean(c_arr))
            median_diffs.append(np.median(t_arr) - np.median(c_arr))
            if ref_kmers is not None:
                contig, pos = key
                contigs.append(contig)
                positions.append(pos)
                ref_kmer_list.append(ref_kmers.kmer_at(contig, pos))
            else:
                ids.append(key)
        if ref_kmers is not None:
            df = pd.DataFrame(data = {'contig' : contigs,
                                      'ref_pos' : positions,
                                      'kmer' : ref_kmer_list,
                                      test_type : test_values,
                                      'pvalue' : p_values,
                                      'mean_diff' : mean_diffs,
                                      'median_diff' : median_diffs,
                                      'n_treated' : n_treated,
                                      'n_control' : n_control})
            sort_cols = ['contig', 'ref_pos']
        else:
            df = pd.DataFrame(data = {id_col : ids,
                                      test_type : test_values,
                                      'pvalue' : p_values,
                                      'mean_diff' : mean_diffs,
                                      'median_diff' : median_diffs,
                                      'n_treated' : n_treated,
                                      'n_control' : n_control})
            sort_cols = [id_col]
        if scaled:
            df['ks_scaled'] = scaled
        df = df.sort_values(by = sort_cols)
        return df


    def score_column(self, df, test):
        """The column to threshold and to find peaks on. For KS the raw statistic is
        replaced by D * sqrt(n1 n2 / (n1 + n2)), the argument of the Kolmogorov limit
        law: monotone in the p-value, aware of coverage, and it cannot underflow.

        cvm and ad are used as they come. The Anderson-Darling statistic scipy returns is
        already standardised by its null mean and standard deviation for the observed
        sample sizes and tie pattern, so it is comparable between k-mers without further
        scaling. The Cramer-von Mises statistic is only comparable asymptotically and is
        not tie corrected, so on the discrete data types (qs, sig_dt) its null grows with
        the imbalance between the two sample sizes; prefer ks or ad there."""
        return 'ks_scaled' if test == 'ks' and 'ks_scaled' in df.columns else test


    def decoy_threshold(self, obs, null):
        """Most permissive cut whose counted false discovery rate stays within the target."""
        obs = np.sort(np.asarray(obs))[::-1]
        null = np.sort(np.asarray(null))[::-1]
        n_null = np.searchsorted(-null, -obs, side = 'right')
        fdr = n_null / np.arange(1, len(obs) + 1)
        ok = np.flatnonzero(fdr <= self.calibration_fdr)
        if not len(ok):
            return np.inf, np.nan
        return obs[ok[-1]], float(fdr[ok[-1]])


    def key_columns(self):
        return ['contig', 'ref_pos'] if self.use_refpos else ['kmer']


    def add_call_columns(self, df, test, null_df = None):
        """Add the FDR decision and, in ref_pos mode, the length of the run of
        neighbouring significant positions that each significant position belongs to."""
        if len(df) == 0 or 'pvalue' not in df.columns:
            return df
        df = df.copy()
        if test == 'ad' and self.verb:
            # scipy interpolates the Anderson-Darling p value from a small table and
            # clamps it at the ends, so it saturates long before the statistic does.
            # The BH decision below is therefore far coarser than for ks or cvm.
            print("Note: Anderson-Darling p values are clamped by scipy, so the",
                  "'significant' column is coarse for ad. Use 'called', which thresholds",
                  "the statistic against the calibration control.")
        df['significant'], _ = fdrcorrection(df['pvalue'], alpha = self.alpha)
        if self.use_refpos and 'ref_pos' in df.columns:
            df = self.add_run_lengths(df)
        # Primary rule: the control-versus-control threshold. Secondary rule: the FDR
        # decision, narrowed by the run length filter when the positions are adjacent.
        if null_df is not None and len(null_df):
            keys = self.key_columns()
            score = self.score_column(df, test)
            df = df.merge(null_df[keys + [score]].rename(columns = {score : f'null_{score}'}),
                          on = keys, how = 'left')
            threshold, achieved = self.decoy_threshold(df[score].dropna(), null_df[score].dropna())
            df['called'] = df[score] >= threshold
            if self.verb:
                if not np.isfinite(threshold):
                    print(f"Calibrated on {self.calibration_ids}: no {score} cut reaches a counted"
                          f" FDR of {self.calibration_fdr}, so nothing is called. The treated sample"
                          f" is indistinguishable from the calibration control at this depth.")
                else:
                    print(f"Calibrated on {self.calibration_ids} at a target FDR of"
                          f" {self.calibration_fdr}: {score} threshold {threshold:.4g},"
                          f" counted FDR {achieved:.2%}, {int(df['called'].sum())} of {len(df)} called.")
        else:
            df['called'] = df['significant']
            if self.min_run_length > 1 and 'run_length' in df.columns:
                df['called'] = df['called'] & (df['run_length'] >= self.min_run_length)
                if self.verb:
                    print(f"No calibration sample given. FDR at {self.alpha} with a run length"
                          f" filter of {self.min_run_length}: {int(df['called'].sum())} of"
                          f" {len(df)} called.")
        df = self.mark_peaks(df, self.score_column(df, test))
        return df


    def mark_peaks(self, df, test):
        """Mark local maxima of the test statistic among the called positions.

        is_peak    : no called position with a higher statistic lies closer than
                     peak_separation bases (so 2 means the immediate neighbours).
        is_primary : the resolved call. With target_bases set, each peak is moved to
                     the nearest position that could carry the modification and
                     duplicates are merged; without it, the strongest peak within
                     peak_merge_distance wins.
        Ties are resolved in favour of the position closer to the contig start.
        """
        if not self.use_refpos or 'ref_pos' not in df.columns or self.peak_separation < 1:
            # positions are only adjacent in ref_pos mode, so there is nothing to localise
            return df
        df['is_peak'] = False
        df['is_primary'] = False
        for contig, g in df[df['called']].groupby('contig', sort = False):
            g = g.sort_values('ref_pos')
            pos, val, idx = g['ref_pos'].to_numpy(), g[test].to_numpy(), g.index.to_numpy()
            peak = np.zeros(len(pos), dtype = bool)
            for i in range(len(pos)):
                lo = np.searchsorted(pos, pos[i] - self.peak_separation + 1, 'left')
                hi = np.searchsorted(pos, pos[i] + self.peak_separation - 1, 'right')
                w = np.arange(lo, hi)
                w = w[w != i]
                peak[i] = not ((val[w] > val[i]) | ((val[w] == val[i]) & (pos[w] < pos[i]))).any()
            pk = np.flatnonzero(peak)
            df.loc[idx[peak], 'is_peak'] = True
            if self.target_bases is not None:
                # move every peak onto the nearest base that could carry the modification
                cand = self.candidate_positions(df, contig)
                best = {}
                for i in pk:
                    near = cand[np.abs(cand - pos[i]) <= self.snap_distance]
                    if not len(near):
                        continue
                    # nearest candidate, ties resolved toward the contig start
                    q = near[np.lexsort((near, np.abs(near - pos[i])))[0]]
                    if q not in best or val[i] > best[q]:
                        best[q] = val[i]
                keys = df.index[(df['contig'] == contig) & df['ref_pos'].isin(list(best))]
                df.loc[keys, 'is_primary'] = True
            else:
                primary = peak.copy()
                for i in pk:
                    near = pk[(np.abs(pos[pk] - pos[i]) <= self.peak_merge_distance) & (pk != i)]
                    if ((val[near] > val[i]) | ((val[near] == val[i]) & (pos[near] < pos[i]))).any():
                        primary[i] = False
                df.loc[idx[primary], 'is_primary'] = True
        if self.verb:
            how = (f"snapped to the nearest {''.join(sorted(self.target_bases))} within"
                   f" {self.snap_distance}") if self.target_bases is not None else \
                  f"merged within {self.peak_merge_distance}"
            print(f"Peaks among the called positions: {int(df['is_peak'].sum())}"
                  f" (separation {self.peak_separation}); {int(df['is_primary'].sum())}"
                  f" calls after resolving ({how}).")
        return df


    def candidate_positions(self, df, contig):
        """Positions whose sensed base could carry the modification. The sensed base sits
        at index h + sensing_offset of the reported k-mer."""
        slot = self.h + self.sensing_offset
        g = df[df['contig'] == contig]
        keep = [isinstance(k, str) and len(k) > slot and k[slot].upper() in self.target_bases
                for k in g['kmer']]
        return np.sort(g['ref_pos'].to_numpy()[np.array(keep, dtype = bool)])


    def add_run_lengths(self, df):
        df['run_length'] = 0
        for contig, g in df[df['significant']].groupby('contig', sort = False):
            g = g.sort_values('ref_pos')
            pos = g['ref_pos'].to_numpy()
            # a gap of up to run_gap missing positions still counts as the same run
            run_id = np.cumsum(np.concatenate(([True], np.diff(pos) > self.run_gap + 1)))
            counts = pd.Series(run_id).value_counts()
            df.loc[g.index, 'run_length'] = pd.Series(run_id).map(counts).to_numpy()
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
        if len(df) == 0 or 'kmer' not in df.columns:
            return
        if 'is_primary' in df.columns and df['is_primary'].any():
            filt = df['is_primary'].to_numpy()
        elif 'called' in df.columns:
            filt = df['called'].to_numpy()
        elif 'significant' in df.columns:
            filt = df['significant'].to_numpy()
        else:
            filt, _ = fdrcorrection(df['pvalue'], alpha = self.alpha)
        df_f = df[filt]
        if self.use_refpos:
            # Name the sequences by their reference position so that the hits can be traced back
            names = [f"{c}:{p}" for c, p in zip(df_f.contig, df_f.ref_pos)]
        else:
            names = None
        write_out_sequences(list(df_f.kmer),
                            os.path.join(self.kmer_out_path, "_".join([self.run_id, self.data_type, test]) + '.fasta'),
                            names = names)


    def compare_samples(self):
        if self.n_splits > 1:
            kmer_set = get_all_kmers(k=self.k)
            split_size  = len(kmer_set) // self.n_splits
        results, nulls = {}, {}
        for split in range(self.n_splits):
            if self.n_splits > 1:
                kmers = kmer_set[split_size*split:split_size*(split+1)]
            else:
                kmers = None
            control_data = self.load_distribution(self.control_ids, kmers)
            treated_data = self.load_distribution(self.treated_ids, kmers)
            calib_data = None
            if self.calibration_ids is not None:
                calib_data = self.load_distribution(self.calibration_ids, kmers)
                if calib_data is None:
                    print(f"Warning: no values found for the calibration sample"
                          f" {self.calibration_ids}. Falling back to the FDR rule.")
            for tt in self.test_types:
                test_df = self.test_samples(control_data, treated_data, tt)
                # The same test between the two controls, where nothing is modified
                # The calibration pair is two controls, whose own scale factor is 1.0000,
                # so it is tested uncorrected. Applying the treated sample's coefficients
                # here would distort the null the threshold is read from.
                null_df = self.test_samples(control_data, calib_data, tt, correct = False) \
                          if calib_data is not None else None
                if split >= 1:
                    results[tt] = pd.concat([results[tt], test_df])
                    if null_df is not None:
                        nulls[tt] = pd.concat([nulls[tt], null_df])
                else:
                    results[tt] = test_df
                    nulls[tt] = null_df
        for test, df in results.items():
            df = self.add_call_columns(df, test, nulls.get(test))
            self.save_test_dataframe(df, test)
            self.write_significant_kmers(df, test)
        if self.verb:
            print(f"Finished testing of {self.control_ids} and {self.treated_ids} using {self.test_types}")

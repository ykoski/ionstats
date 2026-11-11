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
                 ref_pos = False,
                 shift = 12,
                 sensing_offset = 2,
                 n_vals = 65536,
                 n_threads = 1,
                 uc4_table_path = None,
                 n_bins = None,
                 normalization = 'iqr',
                 norm_quantiles = None,
                 norm_const = None,
                 norm_table = None,
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
        self.k = k
        self.h = self.k // 2
        # ref_name is only a contig filter, ref_pos selects the mode. Setting ref_name
        # alone still enables ref_pos mode for backwards compatibility.
        self.ref_names = None
        if self.ref_name is not None:
            self.ref_names = set(self.ref_name.split(","))
        self.use_refpos = ref_pos or (self.ref_names is not None)
        # The motor protein sits approximately 12 bases downstream of the sensing region,
        # so the dwell times of sig_shifted_dt are grouped by the sequence that is read
        # shift bases later. Negative values shift upstream instead. This is a k-mer mode
        # feature: in ref_pos mode the dwell times are already in the right coordinates.
        self.shift = int(shift)
        if self.use_refpos and self.data_type == 'sig_shifted_dt' and self.shift != 0:
            print("Warning: shift is not used in ref_pos mode, where the dwell times are",
                  "already grouped by the position they were measured at.",
                  "sig_shifted_dt is identical to sig_dt here.")
            self.shift = 0
        # The base that dominates the ionic current is not the centre of the k-mer: in the
        # r10.4.1 9-mer model the base at index h + sensing_offset does. In ref_pos mode the
        # values are therefore keyed by the position of that base rather than by the k-mer
        # centre, which also makes forward and reverse reads agree on the same position.
        self.sensing_offset = int(sensing_offset)
        self.n_vals = n_vals
        self.signal_coeff = 10 / (self.n_vals-1)
        # --- per-read re-normalisation -----------------------------------
        # uncalled4 writes uc/ud already z-scored per read (mean 0, std 1).
        # That is not robust: a rare population of deeply blocked events
        # inflates the read's std, so every other event in the read is divided
        # by a larger number and compressed toward the read mean. Pooled over
        # reads this appears as a level-dependent shift -- k-mers above the
        # read mean read lower, k-mers below it read higher -- which is
        # indistinguishable from a genuine pull toward a fixed level.
        #
        # Measured on this dataset: a treated sample's std runs 2.08% above
        # control while its IQR runs only 0.21% above, and the 2.08% figure
        # accounts for the entire apparent effect.
        #
        #   iqr     (median, IQR/1.349)   default; the central 50% cannot be
        #                                 moved by a 0.4% tail population
        #   mad     (median, MAD/0.6745)  Tombo's choice; more robust than std
        #                                 but still responds to the tail
        #   zscore  (mean, std)           reproduces the previous behaviour
        #   none    leave uncalled4's values untouched (identical to zscore
        #           in effect, but skips the work)
        #
        # The divisors are the normal-consistency constants, so for Gaussian
        # data all three agree. Real signal is heavy-tailed, so the ROBUST
        # scales come out larger than the std and the output is on a slightly
        # compressed scale relative to z-scored values. Dimensionless results
        # (regression slopes, R^2) are unaffected; absolute levels and any
        # hard thresholds must be re-derived. Set norm_const to override the
        # divisor if you need to match an existing scale exactly.
        self.normalization = str(normalization).lower()
        if self.normalization not in ('iqr', 'mad', 'zscore', 'none'):
            raise ValueError(f"normalization must be iqr, mad, zscore or none; "
                             f"got {normalization!r}")
        # Which quantile range the 'iqr' method uses. 25-75 is the IQR default
        self.norm_quantiles = tuple(norm_quantiles) if norm_quantiles else (25.0, 75.0)
        if not (0 <= self.norm_quantiles[0] < self.norm_quantiles[1] <= 100):
            raise ValueError(f"norm_quantiles must be 0 <= lo < hi <= 100; "
                             f"got {self.norm_quantiles}")
        # Optional RAW-SIGNAL normalisation. Tombo and uncalled4 take the scale
        # from the raw signal samples; this collector takes it from the event
        # means in uc.
        self.norm_table = None
        self.norm_kind = None
        if norm_table:
            import numpy as _np
            z = _np.load(norm_table, allow_pickle = False)
            # Two table formats, told apart by their keys rather than by a flag:
            #   raw-signal table  -> iqr_raw / iqr_evt, a unit-free RATIO that
            #                        rescales the read's event-mean scale
            #   shrinkage table   -> `scale`, the shrunk scale ITSELF, already in
            #                        the collected units
            if 'scale' in z.files:
                self.norm_kind = 'shrink'
                self.norm_table = dict(zip([str(x) for x in z['read_ids']],
                                           z['scale'].tolist()))
            else:
                self.norm_kind = 'ratio'
                self.norm_table = dict(zip(
                    [str(x) for x in z['read_ids']],
                    (z['iqr_raw'] / z['iqr_evt']).tolist()))
            print(f"[INFO] raw-signal normalisation: {len(self.norm_table):,}"
                  f" reads from {norm_table}")
        if norm_const is not None:
            self.norm_const = float(norm_const)
        elif self.normalization == 'iqr':
            # normal-consistency: for Gaussian data the scaled range equals sigma,
            # so every method agrees and the output scale does not depend on which
            # range was picked
            from scipy.stats import norm as _norm
            lo, hi = self.norm_quantiles
            self.norm_const = float(_norm.ppf(hi / 100.0) - _norm.ppf(lo / 100.0))
        else:
            self.norm_const = {'mad': 0.6745, 'zscore': 1.0, 'none': 1.0}[
                self.normalization]
        # Arguments that default to null in the config are parsed as strings
        self.save_every = int(save_every) if save_every is not None else None
        self.n_bins = int(n_bins) if n_bins is not None else None
        self.binning = False
        if self.n_bins is not None:
            if self.use_refpos:
                print("Warning: n_bins is not supported in ref_pos mode and is ignored.",
                      "Values are collected exactly, keyed by (contig, ref_pos).")
                self.n_bins = None
            else:
                self.binning = True
                if self.normalization not in ('none', 'zscore'):
                    raise ValueError(
                        "n_bins and a robust normalization cannot be combined: "
                        "binning maps the ORIGINAL integer value range onto bin "
                        "indices, and re-normalising changes that range, so the "
                        "bins would no longer mean what they say. Set "
                        "normalization: none, or n_bins: null.")
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
            # Binning is k-mer mode only, so the first axis is the k-mer index
            return np.zeros((4**self.k, self.n_bins), dtype=np.uint32)
        # Both k-mer and ref_pos mode collect the observed values into a dictionary.
        # The keys are k-mers in k-mer mode and (contig, ref_pos) tuples in ref_pos mode.
        return {}


    def read_norm_params(self, read):
        """(centre, scale) for this read, in the raw uc integer units.

        Always derived from the uc tag -- the signal means -- whichever data
        type is being collected, so sig_mean and sig_std are put on one
        consistent per-read scale.

        Returns None when the read is too short to estimate a spread, or when
        the spread comes out non-positive; the caller then leaves the read's
        values untouched rather than dividing by zero.
        """
        if self.normalization in ('none', 'zscore'):
            return None
        if self.norm_table is not None:
            # raw-signal centre/scale, already in pA. The values being scaled are
            # uncalled4's z-scored ones, so a pA centre cannot be subtracted from
            # them -- only the SCALE transfers, as a ratio against this read's own
            # event-mean spread. Reads missing from the table fall through to the
            # event-mean estimate rather than being dropped.
            hit = self.norm_table.get(read.query_name)
            if hit is not None and np.isfinite(hit) and hit > 0:
                if self.norm_kind == 'shrink':
                    # `hit` is the scale. The centre stays per-read: a noisy
                    # centre shifts a whole read uniformly and pools out across
                    # reads, whereas a noisy scale multiplies every k-mer in the
                    # read by a wrong factor, which is the error worth removing.
                    try:
                        uc = np.asarray(read.get_tag('uc'), dtype = np.float64)
                    except KeyError:
                        return None
                    if uc.size < 100:
                        return None
                    return float(np.median(uc * self.signal_coeff)), float(hit)
                ratio = hit
                try:
                    uc = np.asarray(read.get_tag('uc'), dtype = np.float64)
                except KeyError:
                    return None
                if uc.size >= 100:
                    uc = uc * self.signal_coeff
                    centre = float(np.median(uc))
                    lo, hi = np.percentile(uc, self.norm_quantiles)
                    ev = float(hi - lo) / self.norm_const
                    if np.isfinite(ev) and ev > 0:
                        # `ratio` is IQR_raw(pA)/IQR_evt(pA), which is unit free.
                        # Multiplying the event-mean scale by it gives the RAW
                        # scale expressed in the collected units.
                        return centre, ev * ratio
        try:
            uc = np.asarray(read.get_tag('uc'), dtype = np.float64)
        except KeyError:
            return None
        if uc.size < 100:
            return None
        # AMPLITUDE units, matching the values after kmer_idx_to_amps
        uc = uc * self.signal_coeff
        centre = float(np.median(uc))
        if self.normalization == 'iqr':
            lo, hi = np.percentile(uc, self.norm_quantiles)
            scale = float(hi - lo)
        else:                                   # mad
            scale = float(np.median(np.abs(uc - centre)))
        if not np.isfinite(scale) or scale <= 0:
            return None
        return centre, scale / self.norm_const


    def apply_normalization(self, values, params):
        """Re-standardise one read's values.

        sig_mean is a location, so centre and scale both apply.
        sig_std is a spread, so only the scale applies, subtracting the
        centre from a standard deviation would be meaningless.
        sig_dt and qs are counts and are never touched.

        Applied after kmer_idx_to_amps, so both the values and the parameters
        are in amplitude units. AEAD types are deliberately not handled here,
        see get_read_values.
        """
        if params is None:
            return values
        centre, scale = params
        if self.data_type == 'sig_mean':
            return ((np.asarray(values, dtype = np.float64) - centre)
                    / scale).astype(np.float16)
        if self.data_type == 'sig_std':
            return (np.asarray(values, dtype = np.float64)
                    / scale).astype(np.float16)
        return values


    def get_read_values(self, read, read_mapping):
        if self.data_type == 'qs' and self.use_refpos:
            # In ref_pos mode the quality scores are read straight off the aligned pairs,
            # so no k-mer window and no reference sequence is needed
            values, read_mapping = self.get_refpos_quality_scores(read)
        elif self.data_type == 'qs':
            refname = self.alignment.get_reference_name(read.reference_id)
            pairs = read.get_aligned_pairs(matches_only = True)
            al_start = read.reference_start
            al_end = read.reference_end
            # query_qualities is in the same orientation as the aligned query sequence,
            # so it can be indexed with the query positions of the aligned pairs
            qs = read.query_qualities
            contig_len = self.len_dict[refname]
            qs_seq_start = max(0, al_start - self.h)
            qs_seq_end = min(contig_len, al_end + self.h)
            # Fetch new sequence that is always in forward orientation and
            # is adjusted for qs_seq_start and qs_seq_end
            qs_sequence = self.reference.fetch(refname, qs_seq_start, qs_seq_end)
            qs_kmers = get_kmers(qs_sequence, self.k)
            # The first k-mer of qs_kmers is centred h bases into the fetched sequence
            values, read_mapping = self.get_matched_quality_scores(pairs, qs, qs_kmers, qs_seq_start + self.h)
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
        norm = self.read_norm_params(read)
        values = self.kmer_idx_to_amps(values)
        if self.data_type == 'AEAD_mean' or self.data_type == 'AEAD_std':
            if self.use_refpos:
                ref_kmers = self.get_kmers_or_refpos(read, force_kmers=True)
            else:
                ref_kmers = read_mapping
            if self.data_type == 'AEAD_mean':
                exp = get_expected_values_for_read(ref_kmers, self.uc4_means)
            else:
                exp = get_expected_values_for_read(ref_kmers, self.uc4_stds)
            # AEAD cannot be made fully robust here.
            #
            # The uc4 table is in uncalled4's OWN per-read z-score convention,
            # so the difference must be taken before any re-normalisation,
            # an expectation on one scale cannot be subtracted from an
            # observation on another.
            diffs = np.abs(exp - values)
            values = compute_AEAD_values(diffs, k = self.k)
            if norm is not None:
                values = (np.asarray(values, dtype = np.float64)
                          / norm[1]).astype(np.float16)
        else:
            values = self.apply_normalization(values, norm)
        return values, read_mapping


    def get_matched_quality_scores(self, pairs, quals, kmers, first_kmer_refpos, max_qs=None):
        # first_kmer_refpos is the reference position that kmers[0] is centred on
        values = []
        out_kmers = []
        for qp, rp in pairs:
            kmer_idx = rp - first_kmer_refpos
            if 0 <= kmer_idx < len(kmers):
                q = quals[qp]
                if max_qs is not None and q > max_qs:
                    q = max_qs
                values.append(q)
                out_kmers.append(kmers[kmer_idx])
        return values, out_kmers


    def apply_shift(self, values, read_mapping):
        # The values are in read order, so a positive shift attributes the value observed
        # at read position i to the k-mer or reference position at read position i + shift,
        # which is the sequence at the motor protein while position i is being sensed
        if self.shift > 0:
            return values[:len(values)-self.shift], read_mapping[self.shift:]
        elif self.shift < 0:
            return values[-self.shift:], read_mapping[:len(read_mapping)+self.shift]
        return values, read_mapping


    def get_refpos_quality_scores(self, read, max_qs = None):
        # query_qualities is in the same orientation as the aligned query sequence,
        # so it can be indexed directly with the query positions of the aligned pairs
        quals = read.query_qualities
        values, refpos = [], []
        for qp, rp in read.get_aligned_pairs(matches_only = True):
            q = quals[qp]
            if max_qs is not None and q > max_qs:
                q = max_qs
            values.append(q)
            refpos.append(rp)
        return values, refpos


    def get_kmers_or_refpos(self, read, force_kmers = False):
        ur = read.get_tag('ur')
        rev = read.is_reverse
        refname = self.alignment.get_reference_name(read.reference_id)
        if not self.use_refpos or force_kmers:
            sequence = self.reference.fetch(refname, ur[0], ur[1])
            if rev:
                sequence = reverse_complement(sequence)
            kmers = get_kmers(sequence, self.k)
            return kmers
        else:
            refpos = np.arange(ur[0]+self.h, ur[1]-self.h)
            if rev:
                refpos = refpos[::-1]
            if self.sensing_offset:
                # Shift along the read, so that the reverse strand moves the other way
                refpos = refpos - self.sensing_offset if rev else refpos + self.sensing_offset
            return refpos
        
    
    def kmer_idx_to_amps(self, values):
        values = np.asarray(values)
        if self.binning:
            # The signal values are signed, so they are shifted to be non-negative
            values = (values + self.n_vals // 2) // self.bw
            if self.data_type != 'AEAD_mean' and self.data_type != 'AEAD_std':
                # Binned values are used as column indices of the count matrix
                return values.astype(np.uint32)
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
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            refname = self.alignment.get_reference_name(read.reference_id)
            if self.ref_names is not None and refname not in self.ref_names:
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
                read_mapping = read_mapping[self.h:len(read_mapping)-self.h]
            elif self.data_type == 'sig_shifted_dt':
                values, read_mapping = self.apply_shift(values, read_mapping)
            assert len(read_mapping) == len(values), f"Mapping should be the same length as values {len(read_mapping)} != {len(values)}"
            if self.use_refpos:
                # The contig name is part of the key, so positions of different
                # contigs are kept apart
                for i in range(len(values)):
                    add_to_dict(data, (refname, int(read_mapping[i])), values[i])
            elif self.binning:
                for i in range(len(values)):
                    kmer_idx = self.kmer_map[read_mapping[i]]
                    data[kmer_idx, values[i]] += 1
            else:
                for i in range(len(values)):
                    add_to_dict(data, read_mapping[i], values[i])
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
        if self.binning and (self.data_type not in ['AEAD_mean', 'AEAD_std']):
            fname += ".npy"
            np.save(fname, data)
        else:
            trim_dict(data)
            fname += ".pickle"
            with open(fname, 'wb') as handle:
                pickle.dump(data, handle, protocol = pickle.HIGHEST_PROTOCOL)
        if self.verb:
            print(f"Saved collector output to {fname}")
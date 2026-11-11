import numpy as np
import os
import itertools
import errno
import pathlib
import yaml
import pysam


CHUNK_SIZE = 1000

def get_read_ids(lst_file):
    read_ids = []
    with open(lst_file, 'r') as lst:
        for line in lst:
            read_ids.append(line.strip())
    return set(read_ids)


def reverse_complement(seq, comps = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N'}):
    return "".join([comps[s] for s in seq[::-1]])


def reformat_ul(ul):
    ul = list(ul)
    while ul[1] < 0:
        ul[1] += ul[0]
        ul = ul[1:]
    return np.array(ul)


def get_kmer_mean_and_std(uc4_table_path):
    d = np.load(uc4_table_path)
    kmers = d['kmer']
    means = d['current.mean']
    stds = d['current_sd.mean']
    mean_dict = {kmers[i] : means[i] for i in range(len(kmers))}
    std_dict = {kmers[i] : stds[i] for i in range(len(kmers))}
    print(f"len kmers in uncalled4 table: {len(kmers)}")
    print(len(kmers[0]))
    return mean_dict, std_dict


def get_expected_values_for_sequence(seq, mean_dict, std_dict, kmer_len = 9):
    h = kmer_len // 2
    seqlen = len(seq)
    means = np.zeros(seqlen - 2*h)
    stds = np.zeros(seqlen - 2*h)
    for i in range(h,seqlen-h):
        kmer = seq[i-h:i+h+1]
        means[i-h] = mean_dict[kmer]
        stds[i-h] = std_dict[kmer]
    return means, stds


def get_expected_values_for_read(seq, exp_dict):
    # The k-mer table covers ACGT only, so any k-mer overlapping an N in the reference
    # has no expected value. Those positions become NaN rather than raising: the AEAD
    # windows that touch them are then dropped by add_to_dict, which leaves the rest of
    # the read intact. A single N in a reference used to abort the whole sample.
    values = np.full(len(seq), np.nan)
    for i, m in enumerate(seq):
        values[i] = exp_dict.get(m, np.nan)
    return values


def add_to_dict(kmer_dict, kmer, value, dtype = np.float16):
    """Dynamically add values to a NumPy array stored in a dictionary."""
    if kmer is None:
        return
    # NaN reaches here only from an AEAD window that overlapped a k-mer missing from the
    # expected value table. Dropping it keeps the stored distributions clean. The test is
    # the cheap self-comparison rather than np.isfinite, since this runs once per value.
    if value != value:
        return
    if kmer in kmer_dict:
        arr, length = kmer_dict[kmer]
        if length == arr.size:  # Need to expand
            new_size = arr.size + CHUNK_SIZE
            new_arr = np.empty(new_size, dtype=arr.dtype)
            new_arr[:length] = arr  # Copy existing values
            arr = new_arr
        
        arr[length] = value  # Assign new value
        kmer_dict[kmer] = (arr, length + 1)  # Update length
    else:
        arr = np.empty(CHUNK_SIZE, dtype=dtype)
        arr[0] = value
        kmer_dict[kmer] = (arr, 1)  # Start with 1 element


def get_kmers(sequence, k, of = True):
    # of : only full, meaning that the k-mers at the edges
    # of the sequence are skipped
    kmers = []
    halfk = k // 2
    seqlen = len(sequence)
    for i in range(seqlen):
        if of and (i < halfk or i >= seqlen - halfk):
            continue
        if i < halfk:
            kmers.append((halfk-i) * "-" + sequence[:i+halfk+1])
        elif i >= seqlen - halfk:
            kmers.append(sequence[i-halfk:] + (halfk - (seqlen - 1 - i)) * "-")
        else:
            kmers.append(sequence[i-halfk:i+halfk+1])
    return kmers


def get_all_kmers(k, bases = ['A', 'C', 'G', 'T']):
    return [''.join(p) for p in itertools.product(bases, repeat=k)]


def trim_dict(dictionary):
    for kmer in dictionary:
        arr, length = dictionary[kmer]
        dictionary[kmer] = arr[:length]  # Remove unused capacity


def compute_AEAD_values(diffs, k = 9):
    h = k // 2
    aeads = []
    for i in range(h, len(diffs)-h):
        aeads.append(np.mean(diffs[i-h:i+h+1]))
    return aeads


def update_dict(d1, d2, kmer_set = None, np_arrays = True):
    newd = {}
    k1 = set(d1.keys())
    k2 = set(d2.keys())
    if kmer_set is not None:
        k1 = k1 & set(kmer_set)
        k2 = k2 & set(kmer_set)
    for key in k1 & k2:
        if np_arrays:
            newd[key] = np.concatenate((d1[key], d2[key]))
        else:
            newd[key] = d1[key] + d2[key]
    for key in k1 - k2:
        newd[key] = d1[key]
    for key in k2 - k1:
        newd[key] = d2[key]
    return newd


def counts_to_array(counts, conversion_vals):
    array = []
    start, binsize = conversion_vals
    for i, c in enumerate(counts):
        array.extend([start + i*binsize]*int(c))
    return np.array(array)


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
        print(f"Creating directory: {path}")
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print(f"Path {path} already exists.")
            raise


def get_matching_files(data_path, strand, data_type, ids):
    matching_files = []
    all_files = os.listdir(data_path)
    for id in ids:
        for f in all_files:
            full_path = os.path.join(data_path, f)
            if strand != 'both':
                if id in f and strand in f and data_type in f:
                    matching_files.append(full_path)
            else:
                if id in f and data_type in f:
                    matching_files.append(full_path)
    return matching_files


def get_file_path(id, bpath):
    for root, dirs, files in os.walk(bpath):
        for f in files:
            if id in f:
                return os.path.join(root, f)
            

def write_out_sequences(seqs, filename, names = None):
    with open(filename, 'w') as ofile:
        for e, seq in enumerate(seqs):
            if "-" not in seq:
                name = f"kmer:{e}" if names is None else names[e]
                ofile.write(f">{name}\n")
                ofile.write(seq + "\n\n")
    return


def get_id_column(line_list, col, value):
    try:
        id_index = line_list.index(col)
    except ValueError:
        print(f"{value} column {col} is not in the header!",
                "Check the sequencing summary file and 'filter_column' value. Exiting...")
        exit()
    return id_index


def overlay_config(base, new, prefix = ""):
    """Write the values of new over base, recursing into nested sections."""
    for key, val in new.items():
        if isinstance(val, dict) and isinstance(base.get(key), dict):
            overlay_config(base[key], val, f"{prefix}{key}.")
        else:
            base[key] = val


def load_config(config_path = None, default_config_path = "default_config.yaml"):
    # A custom config is layered on top of the default one, so that a config written
    # before a new argument existed still runs, with the default filled in for it
    base_dir = pathlib.Path(__file__).resolve().parent.parent
    with open(os.path.join(base_dir, default_config_path), 'r') as f:
        config = yaml.safe_load(f)
    if config_path is None:
        return config
    with open(pathlib.Path(config_path), 'r') as f:
        custom = yaml.safe_load(f) or {}
    filled = sorted(set(config) - set(custom))
    overlay_config(config, custom)
    if filled:
        print(f"[INFO] Arguments not set in {config_path}, using defaults: {', '.join(filled)}")
    return config
    

def infer_type(val):
    if isinstance(val, bool):
        # Think about store_true/store_false option
        return "bool"
    elif isinstance(val, int):
        return int
    elif isinstance(val, float):
        return float
    elif isinstance(val, list):
        return 'list'
    else:
        return str   
    

def get_reference_kmers(ref_path, ref_name, k):
    reference = pysam.FastaFile(ref_path)
    ref_seq = reference.fetch(reference = ref_name)
    ref_kmers = get_kmers(ref_seq, k, of = False)
    return ref_kmers


class RefKmerFetcher():
    """Fetch the reference k-mer that was in the pore at a (contig, position) pair.

    The position is the base that dominates the signal, which sits at index
    h + sensing_offset of the k-mer rather than at its centre, so the fetched
    window is shifted by -sensing_offset with respect to the position.

    One contig sequence is cached at a time, which is enough when the positions
    are iterated in sorted order, and avoids materialising the k-mers of every
    contig of the reference.
    """

    def __init__(self, ref_path, k, sensing_offset = 0):
        self.reference = pysam.FastaFile(str(ref_path))
        self.k = k
        self.h = k // 2
        self.sensing_offset = sensing_offset
        self.contig = None
        self.sequence = None

    def kmer_at(self, contig, pos):
        if contig != self.contig:
            self.sequence = self.reference.fetch(reference = contig)
            self.contig = contig
        seqlen = len(self.sequence)
        start, end = pos - self.h - self.sensing_offset, pos + self.h - self.sensing_offset + 1
        # Pad the k-mers at the contig edges, as get_kmers does
        left = "-" * max(0, -start)
        right = "-" * max(0, end - seqlen)
        return left + self.sequence[max(0, start):min(end, seqlen)] + right

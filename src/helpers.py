import numpy as np
import os
import itertools
import errno
import pathlib
import yaml


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
    values = np.zeros(len(seq))
    for i, m in enumerate(seq):
        values[i] = exp_dict[m]
    return values


def add_to_dict(kmer_dict, kmer, value, dtype = np.float16):
    """Dynamically add values to a NumPy array stored in a dictionary."""
    if kmer is None:
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


def get_file_path(id, bpath):
    for root, dirs, files in os.walk(bpath):
        for f in files:
            if id in f:
                return os.path.join(root, f)
            

def write_out_sequences(seqs, filename):
    with open(filename, 'w') as ofile:
        for e, seq in enumerate(seqs):
            if "-" not in seq:
                ofile.write(f">kmer:{e}\n")
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


def load_config(config_path = None, default_config_path = "default_config.yaml"):
    base_dir = pathlib.Path(__file__).resolve().parent.parent
    default_config_path = os.path.join(base_dir, default_config_path)
    pth = pathlib.Path(config_path or default_config_path)
    with open(pth, 'r') as f:
        return yaml.safe_load(f)
    

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
    

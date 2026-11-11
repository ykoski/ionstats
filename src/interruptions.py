import pysam
from helpers import *

class IonStatsReadInterruptions():

    def __init__(self,
                 bam_path,
                 summary_path,
                 sample_id,
                 ref_path,
                 out_path,
                 run_id,
                 summary_id = None,
                 read_id_list = None,
                 window = 25,
                 skip_edges = True,
                 use_padding = False,
                 write_ids = True,
                 filter_column = 'barcode_rear_score',
                 read_id_column = 'read_id',
                 id_column = 'barcode_arrangement',
                 filter_threshold = 60,
                 invert_filter = False):
        self.summary_path = summary_path
        self.alignment = pysam.AlignmentFile(bam_path)
        self.reference = pysam.FastaFile(ref_path)
        self.out_path = out_path
        self.sample_id = sample_id
        # The sequencing summary identifies a sample by its barcode, while the bam files
        # are usually named after the sample. summary_id lets the two differ; without it
        # the sample id would have to be the barcode to match the summary and the sample
        # name to find the bam, which cannot both hold.
        self.summary_id = summary_id if summary_id is not None else sample_id
        # A file of read ids to use instead of filtering the summary. This is what lets
        # the pod5 signal filter sit in the middle of the chain: the summary gives the
        # interrupted reads, the filter splits them into blocking and passing by their
        # end signal, and each subset comes back here for sequence extraction.
        self.read_id_list = read_id_list
        self.run_id = run_id
        self.read_ids = []
        self.sequences = []
        self.read_id_index = []
        self.write_ids = write_ids
        self.filter_column = filter_column
        self.read_id_column = read_id_column
        self.id_column = id_column
        self.filter_threshold = filter_threshold
        # Normally the reads kept are those below the threshold, i.e. the interrupted
        # ones. Inverting keeps the rest, which gives the fully sequenced read ends,
        # the natural control set for a motif comparison, since they are the same
        # construction from the same library.
        self.invert_filter = invert_filter
        self.window = window
        # A read whose mapping end sits near a contig boundary has its window clipped, so
        # the extracted sequence comes out short. On a fragmented draft assembly that is
        # common, and the truncated sequences carry no information about why the read
        # ended, they are where the aligner ran out of reference. skip_edges drops them,
        # use_padding keeps them padded with '-' instead. Silently returning a short
        # sequence, which is what this did before, is the one option that hides the issue.
        self.skip_edges = skip_edges
        self.use_padding = use_padding
        self.n_edge = 0
        self.verb = True
    

    def get_interrupted_read_ids(self):
        if self.read_id_list is not None:
            with open(self.read_id_list, 'r') as ifile:
                self.read_ids = [line.strip() for line in ifile if line.strip()]
            self.read_id_set = set(self.read_ids)
            self.read_id_pos = {r: i for i, r in enumerate(self.read_ids)}
            if self.verb:
                print(f"Read {len(self.read_ids)} read ids from {self.read_id_list};"
                      f" the sequencing summary is not used.")
            return
        with open(self.summary_path, 'r') as ifile:
            for i, line in enumerate(ifile):
                # strip first: without it the newline stays on the last field, so a
                # column that happens to be last can never be matched by name
                line_list = line.rstrip("\n").split("\t")
                if i == 0:
                    # Header line
                    id_index = get_id_column(line_list, self.id_column, 'Id')
                    read_id_index = get_id_column(line_list, self.read_id_column, 'Read id')
                    filter_index = get_id_column(line_list, self.filter_column, 'Filter')
                else:
                    if line_list[id_index] == self.summary_id:
                        try:
                            filter_value = float(line_list[filter_index])
                        except ValueError:
                            # a missing or non-numeric score cannot be classified either way
                            continue
                        keep = (filter_value >= self.filter_threshold if self.invert_filter
                                else filter_value < self.filter_threshold)
                        if keep:
                            self.read_ids.append(line_list[read_id_index])
        self.read_id_set = set(self.read_ids)
        self.read_id_pos = {r: i for i, r in enumerate(self.read_ids)}


    def write_read_ids(self, opath):
        with open(opath, 'w') as ofile:
            for read_id in self.read_ids:
                ofile.write(read_id + '\n')


    def write_sequences(self, opath):
        with open(opath, 'w') as ofile:
            for i in range(len(self.sequences)):
                read_id = self.read_ids[self.read_id_index[i]]
                sequence = self.sequences[i]
                ofile.write(">" + read_id + '\n')
                ofile.write(sequence + "\n")
    

    def get_interrupted_sequences(self):
        for read in self.alignment.fetch(until_eof = True):
            if read.is_secondary or read.is_supplementary:
                continue
            rev = read.is_reverse
            if read.query_name in self.read_id_set:
                if rev:
                    end_pos = read.reference_start
                else:
                    end_pos = read.reference_end - 1
                ref_len = self.reference.get_reference_length(read.reference_name)
                want_start, want_end = end_pos - self.window, end_pos + self.window + 1
                start, end = max(0, want_start), min(ref_len, want_end)
                if start != want_start or end != want_end:
                    self.n_edge += 1
                    if self.skip_edges:
                        continue
                seq = self.reference.fetch(read.reference_name, start, end)
                if self.use_padding and not self.skip_edges:
                    seq = "-" * max(0, -want_start) + seq + "-" * max(0, want_end - ref_len)
                if rev:
                    seq = reverse_complement(seq)
                self.sequences.append(seq)
                self.read_id_index.append(self.read_id_pos[read.query_name])


    def read_interruptions(self):
        self.get_interrupted_read_ids()
        if self.write_ids and self.read_id_list is None:
            opath = os.path.join(self.out_path, "_".join([self.sample_id, self.run_id, "interrupted_reads.lst"]))
            self.write_read_ids(opath)
        self.get_interrupted_sequences()
        if self.verb and self.n_edge:
            what = "skipped" if self.skip_edges else \
                   ("padded" if self.use_padding else "left short")
            print(f"{self.n_edge} reads ended within {self.window} bases of a contig"
                  f" boundary and were {what}.")
        sequence_opath = os.path.join(self.out_path, "_".join([self.sample_id, self.run_id, "interruptions.fasta"]))
        self.write_sequences(sequence_opath)
        print(f"Read interruption extraction finished succesfully.")
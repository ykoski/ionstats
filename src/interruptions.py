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
                 window = 25,
                 write_ids = True,
                 filter_column = 'barcode_rear_score',
                 read_id_column = 'read_id',
                 id_column = 'barcode_arrangement',
                 filter_threshold = 60):
        self.summary_path = summary_path
        self.alignment = pysam.AlignmentFile(bam_path)
        self.reference = pysam.FastaFile(ref_path)
        self.out_path = out_path
        self.sample_id = sample_id
        self.run_id = run_id
        self.read_ids = []
        self.sequences = []
        self.read_id_index = []
        self.write_ids = write_ids
        self.filter_column = filter_column
        self.read_id_column = read_id_column
        self.id_column = id_column
        self.filter_threshold = filter_threshold
        self.window = window
    

    def get_interrupted_read_ids(self):
        with open(self.summary_path, 'r') as ifile:
            for i, line in enumerate(ifile):
                line_list = line.split("\t")
                if i == 0:
                    # Header line
                    id_index = get_id_column(line_list, self.id_column, 'Id')
                    read_id_index = get_id_column(line_list, self.read_id_column, 'Read id')
                    filter_index = get_id_column(line_list, self.filter_column, 'Filter')
                else:
                    if line_list[id_index] == self.sample_id:
                        filter_value = float(line_list[filter_index])
                        if filter_value < self.filter_threshold:
                            self.read_ids.append(line_list[read_id_index])


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
            if read.query_name in self.read_ids:
                if rev:
                    end_pos = read.reference_start
                else:
                    end_pos = read.reference_end - 1
                ref_len = self.reference.get_reference_length(read.reference_name)
                start = max(0, end_pos - self.window)
                end   = min(ref_len, end_pos + self.window + 1)

                seq = self.reference.fetch(read.reference_name, start, end)
                if rev:
                    seq = reverse_complement(seq)
                self.sequences.append(seq)
                self.read_id_index.append(self.read_ids.index(read.query_name))


    def read_interruptions(self):
        self.get_interrupted_read_ids()
        if self.write_ids:
            opath = os.path.join(self.out_path, "_".join([self.sample_id, self.run_id, "interrupted_reads.lst"]))
            self.write_read_ids(opath)
        self.get_interrupted_sequences()
        sequence_opath = os.path.join(self.out_path, "_".join([self.sample_id, self.run_id, "interruptions.fasta"]))
        self.write_sequences(sequence_opath)
        print(f"Read interruption extraction finished succesfully.")
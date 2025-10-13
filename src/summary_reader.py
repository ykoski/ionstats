import pysam
from helpers import *

class IonStatsSummaryReader():

    def __init__(self,
                 summary_path,
                 sample_id,
                 out_path,
                 run_id,
                 extract_columns = ['mean_qscore_template', 'sequence_length_template'],
                 id_column = 'barcode_arrangement'):
        self.summary_path = summary_path
        self.out_path = out_path
        self.sample_id = sample_id
        self.run_id = run_id
        self.read_ids = []
        self.sequences = []
        self.read_id_index = []
        self.id_column = id_column
        self.extract_columns = extract_columns

    
    def extract_values(self):
        column_indices = []
        data = {}
        for ec in self.extract_columns:
            data[ec] = []
        with open(self.summary_path, 'r') as ifile:
            for i, line in enumerate(ifile):
                line_list = line.split("\t")
                if i == 0:
                    # Header line
                    for ec in self.extract_columns:
                        column_indices.append(get_id_column(line_list, ec, ec))
                else:
                    for j, ec in enumerate(self.extract_columns):
                        data[ec].append(line_list[column_indices[j]])
        return data


    def read_summary(self):
        data = self.extract_values()
        for key in data.keys():
            arr = np.array(data[key])
            fname = os.path.join(self.out_path, "_".join([self.sample_id, self.run_id, key + ".npy"]))
            np.save(fname, arr)
        return

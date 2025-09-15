import subprocess as sp
from helpers import *

class IonStatsMotifDiscovery():

    def __init__(self,
                 run_id,
                 sequence_dir,
                 out_dir,
                 alphabet,
                 test_type = None,
                 k = 9,
                 overwrite = True,
                 verb = True):
        self.run_id = run_id
        self.sequence_dir = sequence_dir
        self.out_dir = out_dir
        self.alphabet = alphabet
        self.test_type = test_type
        self.k = k
        self.overwrite = overwrite
        self.verb = verb
        self.control_sequences = get_file_path("control", self.sequence_dir)
        if self.control_sequences is None:
            all_kmers = get_all_kmers(9)
            self.control_sequences = os.path.join(self.sequence_dir, "control_sequences.fasta")
            write_out_sequences(all_kmers, self.control_sequences)
        
        
    def motif_discovery(self):
        command = ["streme"]
        sample_path = get_file_path(self.run_id, self.sequence_dir)
        print(sample_path, self.run_id, self.sequence_dir)
        command.extend(["--p", sample_path])
        out_id = self.run_id
        if self.test_type is not None:
            out_id = "_".join([out_id, self.test_type])
        if self.overwrite:
            command.extend(["--oc", os.path.join(self.out_dir, out_id)])
        else:
            command.extend(["--o", os.path.join(self.out_dir, out_id)])
        command.extend(["--n", self.control_sequences])
        command.extend(["--minw", "3"])
        command.extend(["--maxw", str(self.k)])
        command.extend(["--alph", self.alphabet])
        if self.verb:
            print(f"Running STREME with command {command}.")
        sp.run(command)

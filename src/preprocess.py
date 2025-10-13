import subprocess as sp
import os

class IonStatsPreprocess():

    def __init__(self,
                 bam_path,
                 read_path,
                 ref_path,
                 sample_id,
                 run_id,
                 out_path,
                 verb = True):
        self.bam_path = bam_path
        self.read_path = read_path
        self.ref_path = ref_path
        self.sample_id = sample_id
        self.run_id = run_id
        self.out_path = out_path
        self.verb = verb
        if self.verb:
            print(f"Set up preprocessor for sample {self.sample_id}.")


    def signal_alignment(self):
        command = ["uncalled4", "align"]
        command.extend(["--ref", self.ref_path])
        command.extend(["--reads", self.read_path])
        command.extend(["--bam-in", self.bam_path])
        command.extend(["--bam-out", os.path.join(self.out_path, "_".join([self.sample_id, "signal_alignment", self.run_id]) + ".bam")])
        if self.verb:
            print(f"Running uncalled4 with command {command}.")
        sp.run(command)
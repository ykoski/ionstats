# IonStats

IonStats is a Python software package that can be used to collect values of sequencing variables from nanopore sequencing data and compare two different samples.

Currently, IonStats can collect quality score, ionic current mean, ionic current standard deviation, dwell time, motor-protein-shifted dwell time, AEAD mean, and AEAD standard deviation values. The values can be grouped based on either k-mer context (9 by default) or reference position (within one contig).

IonStats also has the feature of analyzing interrupted reads, which are classified based on the rear_barcode_score value in the sequencing_summary.txt file, provided by the MinKNOW software.

## Installation

You can install IonStats by first cloning into the GitHub repository, then setting up the IonStats environment (either via pip or conda, we recommend using conda), and finally adding `ionstats` to PATH with pip install.

    git clone https://github.com/ykoski/ionstats.git

    cd ionstats
    conda env create -f conda_env.yaml
    conda activate ionstats_env
    pip install -e .

You can verify that the installation was successful by running:

    ionstats --help

## Running IonStats

IonStats can be tested with the provided example dataset. You can download this dataset with the `download_example.sh` script:

    ./download_example.sh

To run IonStats, you can use the following command:

    ionstats --config example_config.yaml

Make sure that you're in the ionstats directory while running the command. The default output directory is ./example_output/ but you can specify another directory with `--out_path`.

This will run all five commands of the software consecutively on the example dataset: `preprocess`, `collect`, `test`, `motif_discovery`, and `interrupted_reads`. These subcommands are described in more detail below.

## Commands

### `preprocess`

### `collect`

### `test`

### `motif_discovery`

### `interrupted_reads`

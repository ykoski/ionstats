# IonStats

IonStats is a Python software package that can be used to collect values of sequencing variables from nanopore sequencing data and compare two different samples.

Currently, IonStats can collect quality score, ionic current mean, ionic current standard deviation, dwell time, motor-protein-shifted dwell time, AEAD mean, and AEAD standard deviation values. The values can be grouped based on either k-mer context (9 by default) or reference position (within one contig).

IonStats also has the feature of analyzing interrupted reads, which are classified based on the rear_barcode_score value in the sequencing_summary.txt file, provided by the MinKNOW software.

Check out our preprint **"Compound-specific DNA adduct profiling with nanopore sequencing and IonStats"**, which describes IonStats in detail on bioRxiv: https://www.biorxiv.org/content/10.1101/2025.10.31.685820v1

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

Make sure that you're in the ionstats directory while running the command. The default output directory is `./example_output/`, but you can specify another directory with `--out_path`.

This will run all seven software commands consecutively on the example dataset: `preprocess`, `read_level_values`, `collect`, `stats`, `test`, `motif_discovery`, and `interruptions`. These subcommands are described in more detail below.

## Commands and arguments

The required arguments that the user must specify are: `control`, `treatment`, `bam_path`, `read_path`, `out_path`, and `ref_path`. If you want to collect read-level values and classify interrupted reads for `control` and `treatment` samples, you must also specify the `summary_path` argument. All arguments are read from a YAML configuration file, which is by default the `default_config.yaml`. Users can create custom config files by copying and modifying the default config, or they can overwrite default argument values via the command line, for example, `--argument value`. 

Here are descriptions of the required arguments:
- `control`: Identifier of the control .bam and read files (for the example data: "control"). All reads that belong to the control group must be in one .bam file inside `bam_path`, and the reads should be ordered similarly.
- `treated`: Identifier of the treated .bam and read files (for the example data: "aa").
- `bam_path`: Path to the directory where the aligned .bam files are located (for the example data: example/bam/).
- `read_path`: Path to the directory where the raw sequencing read files are located (for the example data: example/pod5/).
- `out_path`: Path to the location where all output and intermediate files are created (for the example data: example_output/).
- `ref_path`: Path to the .fasta file of the reference genome, where reads are aligned to (for the example data: example/CEN.PK113-7D_Delft_2012_AEHG00000000.fasta).

By default, IonStats will run all seven commands consecutively. If you wish to run only a specific set of commands, you can do this by modifying the `analyses` argument either in a config file or via the command line argument `--analyses`.

Other general arguments to consider:
- `run_id`: An identifier of the IonStats run, used in the output file names (default: IonStats_run).
- `data_types`: List of which data types to collect and compare, by default using qs, sig_mean, sig_std, and sig_dt. Other options include AEAD_mean, AEAD_std, and sig_shifted_dt.
- `k`: K-mer length. IonStats is designed and tested on k = 9, but if the user wants to group the collected values and test values by shorter or longer k-mer context, it is possible with this argument.
- `ref_name`: Contig name, where to focus the analysis. If this argument is specified, the collected values are grouped by reference position instead of k-mer context. Only reads that are aligned to this contig are considered. Useful for short and simple reference genomes, such as plasmids. **NOTE: There is a bug with this approach, so it does not work at the moment! This is under development and should be fixed soon.**
- `n_bins`: Specifying a value for `n_bins` changes the way that data is collected and tested. Instead of storing the values observed for each k-mer in a dictionary, they are instead stored as a count matrix, where the continuous data is binned into `n_bins`. This feature has not been tested, but it should reduce computational costs while sacrificing some precision.

### `preprocess`

Preprocesses the aligned .bam files and .pod5 files with `uncalled4` and outputs the files into `out_path/signal_alignment/`. Arguments that this command requires are `bam_path` to specify the input .bam location, `read_path` to specify the .pod5 or .fast5 location, and `ref_path` to specify the reference genome file (.fasta). Additionally, `run_id` and sample_ids (`control` and `treated`) are used to name the output files.

### `read_level_values`

Reads the sequencing summary file specified by `summary_path` and collects the distributions specified in `extract_columns`, which includes mean quality score, sequence length, and rear barcode score by default. The data is output as numpy arrays into `out_path/read_level_values/`.

### `collect`

Collects the values of `data_types`, groups them by k-mer context or ref_pos for each treatment into `.pickle` or `.npy` (if `ref_name` or `n_bins` is specified) files. These files are output into `out_path/values/`.

### `stats`

Computes statistics for each k-mer or reference position for each `data_type`, which allows simple comparison of treated and control samples. The calculated statistics are mean, standard deviation, median, lowest 1%, highest 1%, lowest 5%, highest 5%, and number of observations. The data is output as .csv files into `out_path/kmer_stats/`.

### `test`

Compares the distributions of values of each `data_type` for each k-mer or reference position using statistical tests specified in `tests` (by default: ks (Kolmogorov-Smirnov) and cvm (Cramér–von Mises criterion)). The test results are output as .csv files into `out_path/tests/`. Additionally, Benjamini-Hochberg FDR-correction is applied with \alpha = `fdr_alpha` (0.05 by default) to find the k-mers that differ significantly after multiple testing correction. These significant k-mers are written as .fasta to `out_path/tests/`. No k-mers are written out if `ref_name` is set.

### `motif_discovery`

Runs STREME on the significant k-mers discovered in the previous step to search for sequence motifs that were enriched in the significant k-mers, and likely associated with treatment effects. A custom STREME alphabet, stored in `resources/custom_alphabet.txt`, is used to consider the k-mers only in forward-strand orientation. The STREME outputs are written into `out_path/streme_outs/`.


### `interruptions`

This command uses the sequencing summary file to classify interrupted reads based on their `barcode_rear_score` values. The interrupted read ids are written out into `out_path/interrupted_reads` as .lst files and the DNA sequences around the mapping end positions (25 bp up- and downstream by default) are written as .fasta files into the same directory.

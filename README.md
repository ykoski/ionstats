# IonStats

IonStats is a python software package that can be used to collect values of sequencing variables from nanopore sequencing data and compare two different samples.

Currently, IonStats can collect quality score, ionic current mean, ionic current standard deviation, dwell time, motor-protein-shifted dwell time, AEAD mean, and AEAD standard deviation values. The values can be grouped by based on either 9-mer context or reference position (within one contig).

IonStats also has the feature of analyzing interrupted reads, which are classified based on the rear_barcode_score value in the sequencing_summary.txt file, provided by the MinKNOW software. This feature will be implemented soon.

## Installation

## Running IonStats

IonStats can be tested with the provided example dataset.

To run IonStats, you can use the following command:



This will run all five commands of the software consecutively on the example dataset: `preprocess`, `collect`, `test`, `motif_discovery`, and `interrupted_reads`. These subcommands are described in more detail below.

## Commands

### `preprocess`

### `collect`

### `test`

### `motif_discovery`

### `interrupted_reads`
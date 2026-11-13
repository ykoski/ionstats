# IonStats

IonStats is a Python software package that can be used to collect values of sequencing variables from nanopore sequencing data and compare two different samples.

Currently, IonStats can collect quality score, ionic current mean, ionic current standard deviation, dwell time, motor-protein-shifted dwell time, AEAD mean, and AEAD standard deviation values. The values can be grouped based on either k-mer context (9 by default) or reference position (across all contigs of the reference, or a subset of them).

IonStats also has the feature of analyzing interrupted reads, which are classified based on the rear_barcode_score value in the sequencing_summary.txt file, provided by the MinKNOW software.

Check out our preprint **"Compound-specific DNA adduct profiling with nanopore sequencing and IonStats"**, which describes IonStats in detail on bioRxiv: https://www.biorxiv.org/content/10.1101/2025.10.31.685820v1

## Installation

You can install IonStats by first cloning into the GitHub repository, then setting up the IonStats environment (either via pip or conda, we recommend using conda), and finally adding `ionstats` to PATH with pip install.

    git clone https://github.com/ykoski/ionstats.git

    cd ionstats
    conda env create -f conda_env.yaml
    conda activate ionstats_env
    pip install uncalled4
    pip install -e .

`uncalled4` is installed after activating the environment, not from
`conda_env.yaml`, and the order matters. It is published only as a source
package -- no wheels for any platform and no conda package -- so pip compiles
its C++ extension during installation.

You can verify that the installation was successful by running:

    ionstats --help

### uncalled4 is Linux-only, and IonStats only needs it for `preprocess`

`uncalled4` states that it requires GCC and is tested only on Linux, and we have
not been able to build it on macOS.

**This only blocks the `preprocess` step.** IonStats runs `uncalled4` as an
external command and never imports it, so everything else works without it:

* remove `preprocess` from `analyses` in your config, and
* start from BAM files that already carry uncalled4's `uc` and `ud` tags.

Run `preprocess` on a Linux machine, or in a container, and copy the resulting
BAMs over. `collect`, `stats`, `test`, `motif_discovery` and `interruptions` are
pure Python and run anywhere.

On a shared system, set `export PYTHONNOUSERSITE=1` before installing, so that
packages in `~/.local` cannot shadow the ones in the environment.

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
- `ref_pos`: Toggle that groups the collected values by reference position instead of by k-mer context. All contigs of the reference are used, and each value is keyed by a `(contig, ref_pos)` pair, so positions of different contigs are kept apart. The position is that of the sensed base, see `sensing_offset`. Useful for short and simple references, such as plasmids or amplicons. Memory use grows with the number of collected values (roughly 2 bytes per observation, i.e. number of positions × coverage), so `save_every` can be used to write the values out in chunks for very high coverage.
- `ref_name`: Contig name, or a comma-separated list of contig names, to restrict the analysis to. Only reads aligned to these contigs are considered. Specifying `ref_name` without `ref_pos` also enables the reference position mode, for backwards compatibility.
- `shift`: Number of bases that the `sig_shifted_dt` dwell times are shifted by (default: 12). The motor protein is approximately 12 bases downstream of the sensing region, so a dwell time measured while a given base is sensed is grouped by the k-mer context that the pore reads `shift` bases later, which is where the motor protein is at the time of the measurement. A negative value shifts upstream instead, and `shift: 0` makes `sig_shifted_dt` identical to `sig_dt`. The first `shift` values of each read are dropped, since their motor protein context is outside the alignment. Note that the output file names do not include the shift, so use `run_id` to keep the results of different shifts apart. This is a k-mer mode feature, where dwell times from many sequence contexts are pooled into k-mer groups. In `ref_pos` mode every dwell time is already grouped by the position it was measured at, so `sig_shifted_dt` is excluded from `data_types` there.
- `sensing_offset`: Position of the base that dominates the ionic current inside the k-mer, counted from the k-mer centre (default: 2). In the r10.4.1 9-mer model the base at index 6 of the 9-mer explains about 63% of the variance of the expected current and the base at index 5 a further 14%, while the centre base explains under 1%. In `ref_pos` mode the collected values are therefore keyed by the position of that base instead of by the k-mer centre, so that reported positions line up with the modified or variant base, and forward and reverse reads agree on the same position. The offset is applied along the read, so it moves in opposite reference directions on the two strands. Set `sensing_offset: 0` to key the values by the k-mer centre instead. This argument does not affect k-mer mode, where the whole k-mer is reported and no single position has to be chosen.
- `normalization`, `norm_quantiles`, `norm_const`: Per-read re-normalisation of the current values, applied as they are collected (default: `iqr`). uncalled4 writes its `uc` and `ud` values already z-scored per read, to mean 0 and standard deviation 1, and that is not robust: a read carrying a rare population of deeply blocked events has an inflated standard deviation, so every other value in that read is divided by a larger number and compressed toward the read mean. Pooled over reads this appears as a level-dependent shift and is easily mistaken for a genuine effect. The options are `iqr` (median and IQR/1.349, the default), `mad` (median and MAD/0.6745, Tombo's choice, more robust than the standard deviation but still responsive to the tail), and `zscore` or `none` to keep uncalled4's own values. The centre and scale are always estimated from the `uc` tag whichever data type is being collected, so that `sig_mean` and `sig_std` end up on one consistent per-read scale: `sig_mean` is a location and has both centre and scale applied, `sig_std` is a spread and is only rescaled, the AEAD types are only rescaled, and `sig_dt` and `qs` are counts and are never touched. `norm_quantiles` sets the percentile range the `iqr` method uses (default `[25, 75]`). `norm_const` is the divisor applied to the robust scale; `null` uses the normal-consistency constant for the chosen method, so all methods agree on Gaussian data. Real signal is heavy-tailed, so the robust scales come out larger than the standard deviation and the collected values sit on a slightly compressed scale relative to z-scored ones. Dimensionless results such as regression slopes and R^2 are unaffected, but absolute levels, hard thresholds and histogram grids have to be re-derived. A robust normalization cannot be combined with `n_bins`.
- `calibration_control`, `calibration_fdr`, `min_run_length`, `peak_separation`, `peak_merge_distance`, `target_bases`: How the called positions or k-mers are decided. See the `test` section below.
- `treated_scale`, `treated_offset`: Map the treated values onto the control's current scale before testing, as `(x - treated_offset) / treated_scale` (default: `1.0` and `0.0`, i.e. no correction). A treated library can sit on a slightly different current scale than the control: fitting `mean_treated = offset + scale * mean_control` across all k-mers gives a scale of about 0.98 for a strongly affected sample, while two control libraries give 1.00001. Left uncorrected, such a difference produces a per-k-mer effect proportional to the k-mer's own signal level, which is largest for the k-mers with the most extreme current and so manufactures an apparent sequence preference. An affine map does not change the shape of a distribution, so the test is left asking whether the two distributions differ beyond a global rescaling. The calibration pair is always tested uncorrected, since two controls share a scale by construction.
- `jitter`, `js`: Add uniform noise drawn from `[-js, js]` to both samples before the test (default: off, with `js` = 0.01). The discrete data types, `sig_dt` and `qs`, contain many ties, which break the continuity assumption that the KS and CvM null distributions rest on; jitter restores it approximately. The noise is random, but the generator is seeded per comparison, so repeated runs of the same comparison agree. `ad` with the midrank variant is tie-corrected and needs no jitter.
- `n_bins`: Specifying a value for `n_bins` changes the way that data is collected and tested in k-mer mode. Instead of storing the values observed for each k-mer in a dictionary, they are instead stored as a count matrix, where the continuous data is binned into `n_bins`. This feature has not been tested, but it should reduce computational costs while sacrificing some precision. It is not used in `ref_pos` mode, where the exact values are collected instead.

### `preprocess`

Preprocesses the aligned .bam files and .pod5 files with `uncalled4` and outputs the files into `out_path/signal_alignment/`. Arguments that this command requires are `bam_path` to specify the input .bam location, `read_path` to specify the .pod5 or .fast5 location, and `ref_path` to specify the reference genome file (.fasta). Additionally, `run_id` and sample_ids (`control` and `treated`) are used to name the output files.

### `read_level_values`

Reads the sequencing summary file specified by `summary_path` and collects the distributions specified in `extract_columns`, which includes mean quality score, sequence length, and rear barcode score by default. The data is output as numpy arrays into `out_path/read_level_values/`.

### `collect`

Collects the values of `data_types`, groups them by k-mer context or reference position for each treatment into `.pickle` (or `.npy`, if `n_bins` is specified) files. These files are output into `out_path/values/`. A pickled file holds a dictionary of the observed values, keyed by the k-mer in k-mer mode and by the `(contig, ref_pos)` pair in `ref_pos` mode:

    {("5mers_rand_ref_adapter_01", 42): array([0.891, 0.874, ...], dtype=float16), ...}

### `stats`

Computes statistics for each k-mer or reference position for each `data_type`, which allows simple comparison of treated and control samples. The calculated statistics are mean, standard deviation, median, lowest 1%, highest 1%, lowest 5%, highest 5%, and number of observations. The data is output as .csv files into `out_path/kmer_stats/`. The rows are identified by the `kmer` column in k-mer mode, and by the `contig` and `ref_pos` columns in `ref_pos` mode, where the reference k-mer that was in the pore at each position is also reported in the `kmer` column. Note that the position refers to the base at index `k // 2 + sensing_offset` of that k-mer rather than to its centre.

### `test`

#### From test scores to modification calls

After the test statistic is computed, three steps turn it into calls. Each step adds a column to the
test output and nothing is filtered out of it, so any step can be re-derived or re-tuned afterwards.
All three are part of the `test` command and run by default; none of them has to be switched on.

**Step 1 applies in both modes. Steps 2 and 3 are `ref_pos` mode only**, since they rely on the tested
keys being adjacent positions along a contig, which k-mers are not. In k-mer mode the `run_length`,
`is_peak` and `is_primary` columns are not produced at all, and `min_run_length`, `peak_separation`,
`peak_merge_distance` and `target_bases` have no effect; the calls are then the `called` column alone.

**1. Threshold - which positions carry evidence.** At high coverage the p-value stops being a useful
decision rule: it answers "are these two libraries distinguishable here", and two separately prepared
libraries differ at a large fraction of positions.

- `calibration_control`: the sample identifier of a **second control**. The same test is
  run between the two controls, where nothing is modified, and its statistic distribution becomes an
  empirical null that is *counted* rather than assumed. The threshold is the most permissive cut at
  which `(null values above it) / (observed values above it)` stays within `calibration_fdr`
  (default `0.01`) - the same target-decoy logic used in proteomics, so the number means the expected
  fraction of wrong calls rather than an abstract quantile. Works in both k-mer and `ref_pos` mode, and
  the null statistic of every key is written out as `null_<score>` so the cut can be changed without
  rerunning. Treat it as a **lower bound**: a comparison of two controls carries no
  modification-induced degradation, which a genuinely modified sample does.

  For KS the score used for thresholding and for peak finding is not the raw statistic but
  `ks_scaled` = `D * sqrt(n1 n2 / (n1 + n2))`, the argument of the Kolmogorov limit law. The KS
  p-value underflows to exactly zero at high coverage - on the validation panel 26% of positions -
  which leaves a large group BH cannot order, while the scaled statistic is monotone in the p-value,
  aware of coverage, and free of underflow. The per-position observation counts are reported as
  `n_treated` and `n_control`.
- Without a second control the FDR decision is used instead, optionally narrowed by `min_run_length`
  (`ref_pos` mode only): a modification perturbs a whole window, so it appears as a run of neighbouring
  significant positions while library differences tend to give isolated ones. The useful cutoff depends
  on how densely modified the target is - a single modification can only affect `k` consecutive
  positions - so this is a fallback, not a substitute for a control.

Columns: `significant`, `called`, `null_<test>` (only with `calibration_control`), and `run_length`
(`ref_pos` mode only).

**2. Localise - turn the smear into positions.** One modification perturbs every k-mer containing it,
so the called positions form smears rather than points. `is_peak` marks the local maxima of the
statistic among the called positions, using `peak_separation` (default 2, meaning only the immediate
neighbours suppress a candidate). Keeping the separation at 2 matters: in CpG data a CGCG motif can
carry two modifications two bases apart, and a larger separation could never report both.

**3. Resolve - one call per modification.** A single modification usually produces more than one
maximum (on the validation panel, only 15-23% of sites gave exactly one), so the peaks need collapsing.

- With `target_bases` set (e.g. `C` for cytosine methylation, `A` for 6mA): each peak is moved to the
  nearest position that could carry the modification, within `peak_separation - 1` bases, and duplicates
  are merged. This corrects the one-base uncertainty of the peak using the reference, and unlike a
  distance rule it cannot collapse two adjacent modifications. Leave it unset for modifications whose
  target base is not known.
- Without it, the strongest peak within `peak_merge_distance` (default 5) wins. Note that this is close
  to equivalent to simply raising `peak_separation`.

Column: `is_primary`. This is what the significant-k-mer FASTA is written from, so motif discovery
receives roughly one sequence per modification rather than one per affected position. In k-mer mode no
`is_primary` column exists, and the FASTA falls back to the `called` column.

Compares the distributions of values of each `data_type` for each k-mer or reference position using statistical tests specified in `tests` (by default: ks (Kolmogorov-Smirnov), cvm (Cramér–von Mises criterion) and ad (k-sample Anderson-Darling, midrank variant)).

Which test to prefer depends on the data type. `ks` is reported as `ks_scaled`, the statistic multiplied by sqrt(n1 n2 / (n1 + n2)); this is the argument of the Kolmogorov limit law, so its null distribution does not depend on coverage and one threshold applies to every k-mer. `ad` needs no such correction: the statistic scipy returns is already standardised by its null mean and standard deviation for the observed sample sizes, and the midrank variant corrects for ties, which makes it the best behaved test on the discrete data types. `cvm` is only comparable across sample sizes asymptotically and is not tie corrected, so on `qs` and `sig_dt` its null value grows sharply when the two samples differ in depth, and a threshold calibrated on a balanced pair will then call almost everything; prefer `ks` or `ad` in that situation. Note also that scipy clamps the Anderson-Darling p-value, so for `ad` the `called` column (calibrated against a second control) is meaningful while `significant` (Benjamini-Hochberg on the p-value) is coarse. The test results are output as .csv files into `out_path/tests/`, identified by the `kmer` column in k-mer mode and by the `contig` and `ref_pos` columns in `ref_pos` mode. Alongside the test statistic and its p-value, the difference of the means and of the medians between the treated and the control sample are reported as `mean_diff` and `median_diff`. These effect sizes are worth using as the primary ranking: at high coverage the two-sample tests have so much power that the p-value saturates, and a comparison of two unmodified replicates can reach significance at a large fraction of all positions. Additionally, Benjamini-Hochberg FDR-correction is applied with \alpha = `fdr_alpha` (0.05 by default) to find the k-mers that differ significantly after multiple testing correction. These significant k-mers are written as .fasta to `out_path/significant_kmers/`. In `ref_pos` mode, the reference k-mer centred on each significant position is written out and named `>contig:ref_pos`, so that the hits can be traced back to the reference.

Note that in `ref_pos` mode the reads of both strands are pooled at each reference position by default, even though the k-mer context of a position differs between the strands. Use the `strand` argument to analyse the strands separately.

### `motif_discovery`

Runs STREME on the significant k-mers discovered in the previous step to search for sequence motifs that were enriched in the significant k-mers, and likely associated with treatment effects. A custom STREME alphabet, stored in `resources/custom_alphabet.txt`, is used to consider the k-mers only in forward-strand orientation. The STREME outputs are written into `out_path/streme_outs/`.


### `interruptions`

This command uses the sequencing summary file to classify interrupted reads based on their `barcode_rear_score` values. The interrupted read ids are written out into `out_path/interrupted_reads` as .lst files and the DNA sequences around the mapping end positions (25 bp up- and downstream by default) are written as .fasta files into the same directory.

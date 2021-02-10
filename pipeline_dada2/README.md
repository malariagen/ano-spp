# Phylogenetic amplicon sequencing with DADA2

Implementation of amplicon sequencing data processing using dada2 v1.10.0 and additional scripts for data pre-processing and QC.

# Usage

## Setup

This pipeline is based on Snakemake, so all dependencies will be installed in the working directory during pipeline run.

Key input file is `tracking/{sample_set}/samples.csv`. This files includes information on source samples, links to IRODS, and additional metadata.

Mandatory columns are:
- `Replicate` values should be unique for each of tag combinations 
- `Source_sample` for which the genotypes should be reconstructed
- `Run` - Sanger NPG run ID, used to find reads file in IRODS
- `Lane` - Sanger NPG lane ID, used to find reads file in IRODS 
- `Tag` - Sanger NPG tag ID, used to find reads file in IRODS. Note that in case of double tagging, only single tag ID is assigned.

Optional columns frequently used in this project are 
- `Species` of mosquito
- `Sex` of mosquito
- `Extraction` protocol or details
- `Conc_ng_ul` DNA concentration submitted to ampseq 
- `Source_96_well` well in the plate submitted to ampseq
- `Final_PCR_mastermix_concentration` - first 3 runs included tests of different PCR conditions
- `i7_tag` sequence
- `i5_tag` sequence
- `Note` summarizing information on the sample and its processing

## Run pipeline

Submission script is written to work on Sanger farm. The script also has pre-configured `{parent_dir}` containing sample-set-level working directories. This is done to re-use conda across sample sets.

Usage:
```
bash submit.sh {sample_set} {snakemake_params}
```
where:
- `{sample_set}` is the name of sample set name as in `tracking/{sample_set}` directory in this repo
- `{snakemake_params}` - any additional parameters passed to snakemake. Examples: 
	- `-n` - perform a dry run

Example:
```
bash submit.sh run3 -p
```

To run elsewhere, place the FASTQ files in your `{work_dir}` as `import/{replicate}_R{1,2}.fastq.gz`, copy your `samples.csv` to your `{work_dir}`, then run the pipeline with
```
snakemake --use-conda -d {work_dir}
```

During pipeline execution:
- whole pipeline and samples csv are copied to working directory `{parent_dir}/{sample_set}`
- per-replicate crams are downloaded from IRODS and are converted to fastq
- targets are de-multiplexed for each replicate with cutadapt
- haplotype sequences are generated from raw reads for each replicate with dada2
- output data is summarized and filtered, and QC plots are generated

# Key outputs

- `{parent_dir}/{sample_set}/import/{replicate}_R{1,2}.fastq.gz` - reads extracted from IRODS
- `{parent_dir}/{sample_set}/logs/` - LSF logs, check in case of problems
- `{parent_dir}/{sample_set}/dada2/output/stats.tsv` - read counts per target per replicate
- `{parent_dir}/{sample_set}/dada2/output/haplotypes.tsv` - haplotype sequences and statistics
- `{parent_dir}/{sample_set}/dada2/qc/*.pdf` - plots for read counts, filtering, allele count and imbalance, see descriptions below


## QC plot interpretation

All plots summarize information on sample-amplicon combinations. Replicates are combined into samples by `summarise_haplotypes.py` script.

`reads_initial.pdf` - log-scale read counts matching primers. Allows to identify samples/targets failing to produce reads. Important exceptions - P1&P2 should have reads only in the infected mosquitoes. Also, some targets do not work in species distant from _An. gambiae_.

`filter_rate.pdf` - percentage of reads left in the final haplotypes. Usage - similar to previous.

`filter_per_amplicon.pdf`, `filter_per_target.pdf` - filtering per amplicon/target divided by analysis stage.

`allele_counts.pdf` - number of genotypes observed per sample per target. 1 or 2 genotypes are expected from a diploid mosquito. Consistently over 2 alleles in sample suggest contamination, consistently over 2 alleles in target - mispriming/duplication?

`allele_imbalance.pdf` - deviations from 50% in heterozygous markers. This plot is not generated if no heterozygous markers exist.

`allele_freq_cov.pdf` - coverage and allele imbalance distribution. Used to trace low-coverage and low-frequency genotypes.

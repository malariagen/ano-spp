# Phylogenetic amplicon sequencing

Implementation of amplicon sequencing data processing using SeekDeep v2.6.4 and additional scripts for data pre-processing and QC.

# Usage

## Setup

This pipeline is based on Snakemake, so all dependencies will be installed in the working directory during pipeline run.

Key input file is `tracking/{sample_set}/samples.csv`. This files includes information on source samples, links to IRODS, and additional metadata.

Mandatory columns are:
- `Replicate` values should be unique for each of tag combinations. Sanger Sample ID can be used as replicate name. This can be found in SequenceScape manifest.
- `Source_sample` for which the genotypes should be reconstructed. Note that multiple replicates per sample dramatically reduces false positive genotype calls by SeekDeep.
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

Submission script is written to work on Sanger farm4/5. The script also has pre-configured `{parent_dir}` containing all working directories. This is done to re-use conda and singularity across sample sets.

Usage:
```
bash submit.sh {sample_set} {snakemake_params}
```
where:
- `{sample_set}` is the name of sample set name as in `tracking/{sample_set}` directory in this repo
- `{snakemake_params}` - any additional parameters passed to snakemake. Examples: 
	- `-n` - perfrom a dry run
	- `basic_qc` - generate qc results without compressing the analysis folder

Example:
```
bash submit.sh run3 -p
```

To run elsewhere, place the FASTQ files in your `{work_dir}` as `import/{replicate}_R{1,2}.fastq.gz`, copy your `samples.csv` to your `{work_dir}`, then run the pipeline with
```
snakemake --use-conda --use-singularity -d {work_dir}
```

During pipeline execution:
- whole pipeline and samples csv are copied to working directory `{parent_dir}/{sample_set}`
- crams are downloaded from IRODS and are converted to fastq (`import/{sample}_R{1,2}.fastq.gz`)
- SeekDeep is set up and run to generate haplotype sequences from raw reads
- output data is summarized
- SeekDeep analysis directory is compressed
- data QC is performed

# Outputs

- `{parent_dir}/{sample_set}/import/{replicate}_R{1,2}.fastq.gz` - reads extracted from IRODS
- `{parent_dir}/{sample_set}/logs/` - LSF logs, check in case of problems
- `{parent_dir}/{sample_set}/pipeline_seekdeep/` - copy of this pipeline
- `{parent_dir}/{sample_set}/samples.csv` - copy of samples table from tracking
- `{parent_dir}/{sample_set}/seekdeep/analysis.tar.gz` - compressed SeekDeep analysis folder
- `{parent_dir}/{sample_set}/seekdeep/output/extraction.csv` - read counts per target per replicate
- `{parent_dir}/{sample_set}/seekdeep/output/popclustering.csv` - recovered sequences and statistics
- `{parent_dir}/{sample_set}/seekdeep/qc/*.pdf` - plots for read counts, filtering, allele count and imbalance, see descriptions below
- `{parent_dir}/{sample_set}/seekdeep/qc/summary.txt` - lists of failed samples. This file is also copied to `{vector_ampseq_spp_dir}/tracking/{sample_set}/qc_summary.txt`.

## QC plot interpretation

All plots summarize information on sample-amplicon combinations. Replicates are summarized by SeekDeep.

`reads_initial.pdf` - log-scale read counts matching primers. Allows to identify samples/targets failing to produce reads. Important exceptions - P1&P2 should have reads only in the infected mosquitoes. Also, some targets do not work in species distant from An.gambiae.

`reads_final.pdf` - log-scale read counts after filtering. Most importantly, sample-target combinations highlighted in red lost all reads due to filtering.

`filter_rate.pdf` - final reads divided by initial reads. Usage - similar to previous.

`filter_per_amplicon.pdf`, `filter_per_target.pdf` - filtering per amplicon/target divided by SeekDeep stage.

`allele_counts.pdf` - number of genotypes observed per sample per target. 1 or 2 genotypes are expected from a diploid mosquito. Consistently over 2 alleles in sample suggest contamination, consistently over 2 alleles in target - mispriming/duplication?

`allele_imbalance.pdf` - deviations from 50% in heterozygous markers. This plot is not generated if no heterozygous markers exist.

`allele_freq_cov.pdf` - coverage and allele imbalance distribution. Used to trace low-coverage and low-frequency genotypes.

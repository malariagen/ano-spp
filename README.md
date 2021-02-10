# ano-spp
A targeted amplicon sequencing panel to simultaneously identify mosquito species and Plasmodium presence across the entire Anopheles genus 

## Repository structure

- `data` - annotation tables for 62 mosquito and 2 Plasmodium amplicons  
- `pipeline_dada2` - sequencing data processing and QC pipeline based on [DADA2](https://benjjneb.github.io/dada2/tutorial.html)
- `pipeline_seekdeep` - sequencing data processing and QC pipeline based on [SeekDeep](https://seekdeep.brown.edu/)
- `tracking` - sample manifests for pipelines grouped by Illumina MiSeq run with short descriptions of sample sets
- `work` - analyses code and data files; each step has an associated conda environment file `env.yml` listing the dependencies

## Analysis outline

- [1_panel_design](work/1_panel_design) - search for potential amplicon sites in the 21 _Anopheles_ genomes alignment and annotation of the final amplicon set
- [2_plasmodium_rebalancing](work/2_plasmodium_rebalancing) - search for _Plasmodium_ primer concentrations optimal for parasite detection
- [3_plasmodium_qpcr](work/3_plasmodium_qpcr) - comparison of amplicon sequencing and qPCR for _Plasmodium_ detection
- [4_ref_extraction](work/4_ref_extraction) - amplicon sequence extraction from reference genomes (supplementary step)
- [5_synteny_plot](work/5_synteny_plot) - plotting amplicon positions in three mosquito species genomes
- [6_ag1k_extraction](work/6_ag1k_extraction) - amplicon sequence extraction from Ag1000g Phase 2 haplotypes and within-species distances estimation (supplementary step)
- [7_species_id](work/7_species_id) - multiple species dataset exploration, clustering-based species ID, species tree
- [8_ag1k_analysis](work/8_ag1k_analysis) - Ag1000g population structure and diversity based on amplicon sequences
- [9_coi_its](work/9_coi_its) - COI and ITS2 Sanger sequencing data analysis for species ID confirmation, within-species diversity estimates
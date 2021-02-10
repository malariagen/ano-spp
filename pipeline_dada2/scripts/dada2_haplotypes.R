# Perform haplotype inference using dada2 denoising and clustering.
#
# The processing is done for single amplicon, 
# all samples from the run are processed as a batch.

library(dada2)
library(ggplot2)
cat('dada2 version:', as.character(packageVersion("dada2")), '\n', file = stderr())

# multithreading switch 
# set to F to avoid Rcppparallel errors
# set to F to avoid cryptic error in filterAndTrim:
# `'names' attribute [473] must be the same length as the vector [465]`
multithread=FALSE

# status logging function
status <- function(x) cat("##",x,"##\n", file = stderr())

status("Data prep")
# link to snakemake
target <- snakemake@params[["target"]]
rep_names <- snakemake@params[["reps"]]
fq1 <- snakemake@params[["fq1"]]
fq2 <- snakemake@params[["fq2"]]
# plt_q <- snakemake@output[["plt_q"]]
flt_fq1 <- snakemake@params[["flt_fq1"]]
flt_fq2 <- snakemake@params[["flt_fq2"]]
plt_e1 <- snakemake@output[["plt_e1"]]
plt_e2 <- snakemake@output[["plt_e2"]]
haplotypes <- snakemake@output[["haplotypes"]]
statistics <- snakemake@output[["stats"]]

# dataframe of per-sample filenames
# from here on, quality plotting code is commented out 
# d <- data.frame(fq1, fq2, plt_q, flt_fq1, flt_fq2, stringsAsFactors=FALSE)
d <- data.frame(fq1, fq2, flt_fq1, flt_fq2, stringsAsFactors=FALSE)

# exclude samples for which input files don't exist
rep_names <- subset(rep_names, file.exists(d$fq1))
d <- subset(d, file.exists(d$fq1))

if (length(rep_names) == 0) {
	status(paste("No reads originally in amplicon", target))
	file.create(haplotypes)
	file.create(statistics)
	file.create(plt_e1)
	file.create(plt_e2)
} else {

	# print("# plot quality profiles for each sample")
	# for (i in 1:nrow(d)) {
	# 	g = plotQualityProfile(d[i, c("fq1","fq2")])
	# 	ggsave(d[i, "plt_q"], plot=g, width=6, height=3)
	# }

	status("Filter and trim")
	# omit truncLen=c(240,160) - do not remove low quality tail 
	flt_stats = filterAndTrim(
			d$fq1, d$flt_fq1,
			d$fq2, d$flt_fq2,
			maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
			compress=TRUE, multithread=multithread)

	status("Remove samples without reads after filtering")
	d <- subset(d, flt_stats[,"reads.out"] > 0)
	rep_names <- subset(rep_names, flt_stats[,"reads.out"] > 0)

	if (length(rep_names) == 0) {
		# information on samples with all reads filtered out is not saved
		# this will be picked up by summarise_haplotypes
		status(paste("No reads retained after filtering in amplicon", target))
		file.create(haplotypes)
		file.create(statistics)
		file.create(plt_e1)
		file.create(plt_e2) } else {

		status("Compute and plot error profiles")
		# for R1 and R2 separately
		error_profiles1 = learnErrors(d$flt_fq1, multithread=multithread)
		g = plotErrors(error_profiles1)
		ggsave(plt_e1, plot=g)
		error_profiles2 = learnErrors(d$flt_fq2, multithread=multithread)
		g = plotErrors(error_profiles2)
		ggsave(plt_e2, plot=g)

		status("Infer haplotypes")
		# for R1 and R2 separately
		haplotypes1 = derepFastq(d$flt_fq1)
		haplotypes2 = derepFastq(d$flt_fq2)

		status("Infer true haplotypes using denoising-clustering")
		true_haplotypes1 = dada(haplotypes1,
							    err=error_profiles1,
							    verbose=FALSE,
							    multithread=multithread)
		true_haplotypes2 = dada(haplotypes2,
							    err=error_profiles2,
							    verbose=FALSE,
							    multithread=multithread)

		status("Merge R1 and R2 haplotypes")
		merged_true = mergePairs(true_haplotypes1, haplotypes1,
		                         true_haplotypes2, haplotypes2)

		status("Init haplotype matrices") 
		haplotype_matrices = makeSequenceTable(merged_true)

		status("Filter any chimeras")
		# chimeras between samples
		haplotype_matrices.bimera_filtered = removeBimeraDenovo(haplotype_matrices,
			        method="consensus", 
			        multithread=multithread)

		status("Write haplotypes")
		rownames(haplotype_matrices.bimera_filtered) <- rep_names
		write.table(haplotype_matrices.bimera_filtered,
					file=haplotypes,
					col.names=NA, row.names=TRUE,
					sep="\t", quote=FALSE)

		status("Read counts") 
		# convenience function for read counting
		num_unique <- function(x) sum(getUniques(x))

		# single sample output is dada-class
		# multiple sample output is list of dada-classes
		if (length(rep_names) == 1) {
			# do not include samples where all reads were filtered out	
			track <- cbind(subset(flt_stats, flt_stats[,"reads.out"] > 0), 
						num_unique(true_haplotypes1),
						num_unique(true_haplotypes2), 
						num_unique(merged_true),
						rowSums(haplotype_matrices.bimera_filtered))
			} else {
			track <- cbind(subset(flt_stats, flt_stats[,"reads.out"] > 0), 
						sapply(true_haplotypes1, num_unique),
						sapply(true_haplotypes2, num_unique), 
						sapply(merged_true, num_unique),
						rowSums(haplotype_matrices.bimera_filtered))	
			}
		colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
		rownames(track) <- rep_names
		write.table(track,
					file=statistics,
					col.names=NA, row.names=TRUE,
					sep="\t", quote=FALSE)
	}
}
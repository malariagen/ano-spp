'''
- Combine haplotypes and stats across replicates
- Filter sample haplotypes
- QC plots
'''
import itertools
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import gzip
import os
import re

# input files
SAMPLE_FILE = snakemake.input.samples
TARGETS_FILE = snakemake.input.targets
STATS_FILES = snakemake.input.stats
HAPS_FILE = snakemake.input.haplotypes
QC_DIR = snakemake.params.qc_dir
INIT_FREADS_PATTERN = snakemake.params.init_freads_pattern
# output files
STATS_OUT = snakemake.output.stats
HAPS_OUT = snakemake.output.haplotypes
# params
MIN_FRAC = snakemake.params.min_frac
MIN_READS = snakemake.params.min_reads

def msg(msg):
    with open(os.path.join(QC_DIR, 'summary.txt'), 
              mode='a') as o:
        o.write(msg + '\n')

###############
### COMBINE ###
###############

# sample metadata
sample_meta = pd.read_csv(SAMPLE_FILE, dtype='str')

# replicates
all_reps = sample_meta.Replicate.sort_values()

# samples
all_samples = sample_meta.Source_sample.sort_values().unique()

# targets - these can be passed as a list
all_targets = pd.read_csv(TARGETS_FILE, sep='\t', dtype='str')['target'].values

# read counts per replicate
stats_list=list()

def input_read_count(x):
    '''
    Get missing input read counts
    for files with no reads left after filtering
    '''
    # do not re-count positive read counts
    if x.input > 0:
        return x.input
    # file to check - undeclared use of `target`
    fn = INIT_FREADS_PATTERN.format(replicate=x.name, 
                                    target=target)
    # find forward reads file
    if os.path.isfile(fn):
        with gzip.open(fn, 'rb') as f:
            for i, l in enumerate(f):
                pass
        # add a newline at end
        i+=1
        # assert it is a valid fastq
        assert i%4 == 0
        return i//4
    # if no file, return zero reads
    return 0

# get stats for all targets
for sf in STATS_FILES:
    
    # extract target from filename
    target = sf.split('_')[-1].split('.')[0]
    try:
        # read table
        s = pd.read_csv(sf, sep='\t', index_col=0)
        # add missing replicates
        s = s.reindex(all_reps, fill_value=0)
    except:
        # for missing target stats create empty df 
        s = pd.DataFrame(index=all_reps, 
            columns=['input','filtered','denoisedF',
                     'denoisedR','merged','nonchim']).fillna(0)
    s = s.reset_index()
    # add input read counts
    s['input'] = s.apply(input_read_count, axis=1)
    # add target
    s['target'] = target
    stats_list.append(s)
    
    
stats = pd.concat(stats_list)
# manually reset column name
stats.rename(columns={'Replicate':'replicate'}, inplace=True)

# add targets for which no sequencing data was obtained at all
# stats = stats.set_index(['replicate', 'target']) \
#              .reindex(itertools.product(all_reps, all_targets), fill_value=0) \
#              .reset_index()

# haplotypes per replicate
haps_list=list()

# get haplotypes for all targets
for hf in HAPS_FILE:
    # extract target from filename
    target = hf.split('_')[-1].split('.')[0]
    try:
        # read haplotype matrix
        h = pd.read_csv(hf, sep='\t', index_col=0)
        # transform into one hap-rep per line
        h = h.unstack().reset_index()
        # remove zeroes - slow
        h = h[h[0] != 0]
        # add target
        h['target'] = target
        haps_list.append(h)
    # skip empty files generated where no reads for target exist
    except:
        pass
haps = pd.concat(haps_list)
# manually reset column names
haps.columns = ['consensus', 'replicate', 'reads', 'target']

# rep-sample mapping
rep_sample = sample_meta.set_index('Replicate')['Source_sample'].to_dict()
# add to stats
stats['s_Sample'] = stats.replicate.replace(rep_sample)
# add to haps
haps['s_Sample'] = haps.replicate.replace(rep_sample)

# insert per-replicate filtering here. Idea:
# - sequence specific to one replicate

# merge replicate haplotypes
hap_per_sample = haps.groupby(by=['s_Sample', 'target', 'consensus'])['reads'].sum().reset_index()

# add total reads per sample-target combination and fraction of individual haplotype
hap_per_sample['total_reads'] = hap_per_sample.groupby(by=['s_Sample', 'target'])['reads'].transform('sum')
hap_per_sample['frac_reads'] = hap_per_sample['reads'] / hap_per_sample['total_reads']

##############
### FILTER ###
##############

# sample-level filtering: minimal number of reads, minimal fraction cutoffs
hap_per_sample = hap_per_sample[(hap_per_sample.reads >= MIN_READS) &
                                (hap_per_sample.frac_reads >= MIN_FRAC)]

# refresh total and frac read counts
hap_per_sample['total_reads'] = hap_per_sample.groupby(by=['s_Sample', 'target'])['reads'].transform('sum')
hap_per_sample['frac_reads'] = hap_per_sample['reads'] / hap_per_sample['total_reads']

# summarise read counts per sample
read_per_sample = stats.groupby(by=['s_Sample', 'target']).sum()
# add final read counts for filtered haplotypes
read_per_sample['final'] = hap_per_sample.groupby(by=['s_Sample', 'target'])['reads'].sum()
read_per_sample['final'] = read_per_sample.final.fillna(0).astype(int)
read_per_sample = read_per_sample.reset_index()

# write haps
hap_per_sample.to_csv(HAPS_OUT, sep='\t', index=False)
# write stats
read_per_sample.to_csv(STATS_OUT, sep='\t', index=False)

##########
### QC ###
##########

os.makedirs(QC_DIR, exist_ok=True)

def fill_sample_target(df, all_samples=all_samples, all_targets=all_targets):
    '''
    Add NA-only rows and columns to the pandas.DataFrame where
    samples are columns and targets are rows
    '''

    df.columns = df.columns.map(str)
    df.index = df.index.map(str)
    for col in all_samples:
        if col not in df.columns:
            df[col] = np.nan
    for row in all_targets:
        if row not in df.index:
            df = df.append(pd.Series(name=row))
    df = df.sort_index(axis=0)
    df = df.sort_index(axis=1)
    
    return df

def plot_heatmap(data, title, fname, **kwargs):
    if 'cmap' not in kwargs:
        kwargs['cmap']='coolwarm_r'
    
    data = fill_sample_target(data)
    fig_width = data.shape[1] / 4 # samples
    fig_height = data.shape[0] / 5 # targets
    # setting of cbar ax is incompatible with tight layout
    # grid_kws = {'width_ratios': (0.9, 0.02), 'wspace': 0.05}
    # fig, (ax, cbar_ax) = plt.subplots(1, 2, 
    #                                   gridspec_kw=grid_kws, 
    #                                   figsize=(fig_width, fig_height))
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    sns.heatmap(data, 
                ax=ax,
                # cbar_ax=cbar_ax,
                **kwargs)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(QC_DIR, fname), dpi=150)

# initial merged read counts (i.e., matching both primers)
init_reads = read_per_sample \
    .pivot(index='target', columns='s_Sample', values='input')
init_reads = init_reads.replace(0, np.nan)
plot_heatmap(np.log10(init_reads),
             title="Initial reads per sample per target (log10)",
             fname="reads_initial.pdf")

# sucess rate indicates percentages of losses during all stages of extraction and clustering
read_per_sample['filter_rate'] = read_per_sample.final / read_per_sample.input
success_rate = read_per_sample \
    .pivot(index='target', columns='s_Sample', values='filter_rate')
plot_heatmap(success_rate,
             title="Proportion of reads passing all filters per sample per target",
             fname="filter_rate.pdf", 
             center=0.5)

def calc_props(df):
    # calculate proportions removed at every stage
    # note that it is done on the original DataFrame

    df['0_prefilter_prop'] = (df.input - df.filtered) / \
                                        df.input\
    # reverse reads are worse, so report denoising on those
    df['1_noise_prop'] = (df.filtered - df.denoisedR) / \
                                        df.input
    df['2_unmerged_prop'] = (df.denoisedR - df.merged) / \
                                        df.input
    df['3_chimeric_prop'] = (df.merged - df.nonchim) / \
                                        df.input
    df['4_postfilter_prop'] = (df.nonchim - df.final) / \
                                        df.input
    df['5_final_prop'] = df.final / df.input
    return df.loc[:, df.columns.str.endswith('_prop')] \
             .sort_index(axis=1, ascending=False)

# per-amplicon failure reasons
ampl_data = read_per_sample.groupby(by='target').sum()
width = ampl_data.shape[0] / 3
fig, ax = plt.subplots(1, 1, figsize=(12, 3))
frac_ampl = calc_props(ampl_data)
frac_ampl.plot(kind='bar',
               ax=ax,
               cmap='Paired',
               stacked=True,
               title='Read status breakdown per amplicon')
ax.legend([x.split('_')[1] for x in frac_ampl.columns])
plt.savefig(os.path.join(QC_DIR, "filter_per_amplicon.pdf"), dpi=150)

# per-sample failure reasons
sample_data = read_per_sample.groupby(by='s_Sample').sum()
frac_sample = calc_props(sample_data)
split_idx = frac_sample.shape[0] // 2
d1 = frac_sample.iloc[:split_idx, :]
d2 = frac_sample.iloc[split_idx:, :]
width = split_idx / 6
# plot
fig, axs = plt.subplots(2, 1, figsize=(width, 8))
for (i,d) in enumerate([d1, d2]):
    d.plot(kind='bar',
           ax=axs[i],
           cmap='Paired',
           legend=False,
           stacked=True, 
           title=('Read status breakdown per sample' if i==0 else ''))
    # rotate tick labels
#     for tick in axs[i].get_xticklabels():
#         tick.set_rotation(45)
axs[0].legend([x.split('_')[1] for x in frac_sample.columns]);
# fit xticklabels
plt.tight_layout()
plt.savefig(os.path.join(QC_DIR, "filter_per_sample.pdf"), dpi=150)

# allele counts (TODO - discrete color bar)
allele_counts = hap_per_sample.groupby(['s_Sample', 'target'], as_index=False).count()\
    .pivot(index='target', columns='s_Sample', values='consensus') 
max_allele_count = allele_counts.max().max().astype(int) + 1
plot_heatmap(allele_counts,
             title="Alleles per sample per target",
             fname="allele_counts.pdf",
             cmap="coolwarm",
             center=2,
             cbar_kws=dict(ticks=range(max_allele_count)))

# major allele frequency
major_hap_freq = hap_per_sample.groupby(['s_Sample', 'target'], as_index=False)['frac_reads'].max()\
    .pivot(index='target', columns='s_Sample', values='frac_reads')
# ignore perfectly homozygous samples
major_hap_freq = major_hap_freq.replace(1, np.nan)
# plot only if any heterozygous sites exist
if ~major_hap_freq.isna().all(axis=None):
    plot_heatmap(fill_sample_target(major_hap_freq), 
                 title="Major allele imbalance: red - over 0.5, blue - below 0.5",
                 fname="allele_imbalance.pdf",
                 cmap="coolwarm",
                 center=0.5)

# allele frequencies versus read counts
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.scatter(x=hap_per_sample.frac_reads, 
           y=np.log10(hap_per_sample.reads), 
           alpha=0.1)
plt.xlabel('Allele fraction')
plt.ylabel('Read count, log10')
plt.savefig(os.path.join(QC_DIR, "allele_freq_cov.pdf"), dpi=150)


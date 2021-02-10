import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os
import re

SAMPLE_FILE = snakemake.input[0]
EXTRACTION_PROFILE = snakemake.input[1]
POPCLUSTERING_FILE = snakemake.input[2]
TARGETS_FILE = snakemake.input[3]
OUT_DIR = snakemake.params.dir

# threshold for reporting samples with low number of reads:
# proportion of mean totalMatching reads
PROP_READ_COUNT = 0.1
# threshold for reporting samples with 
# low proportion of reads retained after filtering
PROP_RETAINED = 0.1

def plot_heatmap(data, title, fname, **kwargs):
    if 'cmap' not in kwargs:
        kwargs['cmap']='coolwarm_r'
    
#     grid_kws = {'width_ratios': (0.9, 0.03), 'wspace': 0.18}
#     fig, (ax, cbar_ax) = plt.subplots(1, 2, gridspec_kw=grid_kws, figsize=(18,10))
    fig_width = data.shape[1] / 4 # samples
    fig_height = data.shape[0] / 5 # targets
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    sns.heatmap(data, 
                ax=ax, 
                **kwargs)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, fname), dpi=150)

def msg(msg):
    with open(os.path.join(OUT_DIR, 'summary.txt'), 
              mode='a') as o:
        o.write(msg + '\n')

# sample metadata
sample_meta = pd.read_csv(SAMPLE_FILE, dtype={'Source_sample':'str'})
# SeekDeep extraction profile
extr_data = pd.read_csv(EXTRACTION_PROFILE)
# merged SeekDeep popClustering tables
pop_data = pd.read_csv(POPCLUSTERING_FILE)

# read counts prep
read_per_sample = extr_data.groupby(by=['s_Sample', 'target']).sum()
read_per_sample['final'] = pop_data.groupby(by=['s_Sample', 'target'])['c_ReadCnt'].sum()
read_per_sample['final'] = read_per_sample['final'].fillna(0).astype(int)
read_per_sample['failedClustering'] = read_per_sample.good - read_per_sample.final
read_per_sample.drop(columns=['good', 'bad'], inplace=True)

# all samples and targets in the experiment
all_samples = sample_meta['Source_sample'].unique()
all_targets = pd.read_csv(TARGETS_FILE, sep='\t', dtype='str')['target'].unique()

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

# initial merged read counts (i.e., matching at least one primer)
total_reads = read_per_sample \
    .reset_index() \
    .pivot(index='target', columns='s_Sample', values='totalMatching')
plot_heatmap(np.log10(fill_sample_target(total_reads)),
             title="Initial reads per sample per target (log10)",
             fname="reads_initial.pdf")

# sucess rate indicates percentages of losses during all stages of extraction and clustering
success_rate = (read_per_sample.final / read_per_sample.totalMatching) \
    .reset_index() \
    .pivot(index='target', columns='s_Sample', values=0)
plot_heatmap(fill_sample_target(success_rate),
             title="Proportion of reads passing all filters per sample per target",
             fname="filter_rate.pdf", 
             center=0.5)

# final read counts per sample per target
final_reads = read_per_sample.reset_index() \
                             .pivot(index='target', columns='s_Sample', values='final')
# replace zeroes with small number for logscale conversion
final_reads = final_reads.replace(0, 0.009)
# resulting colours:
# - NaN - no reads initially, white
# - 0.009 - all reads removed, red
# - >=1 - some reads retained, grey to blue
plot_heatmap(np.log10(fill_sample_target(final_reads)),
             title="Final reads per sample per target (log10); red - all reads removed",
             fname="reads_final.pdf",
             center=0)

# per-amplicon failure reasons
ampl_data = read_per_sample.groupby(by='target').sum()
for col in ampl_data.drop(columns='totalMatching'):
    ampl_data[col+'_pc'] = ampl_data[col]/ampl_data['totalMatching']
fig, ax = plt.subplots(1, 1, figsize=(18, 3))
ampl_data[['final_pc',
           'failedClustering_pc',
           'failedQuality_pc',
           'failedPairProcessing_pc',
           'failedMinLen_pc',
           'failedMaxLen_pc',
           'failedNs_pc',
           'failedPossibleContamination_pc']] \
                .plot(kind='bar', 
                      stacked=True, 
                      ax=ax, 
                      title='Read status breakdown per amplicon')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "filter_per_amplicon.pdf"), dpi=150)

# per-sample failure reasons
sample_data = read_per_sample.groupby(by='s_Sample').sum()
for col in sample_data.drop(columns='totalMatching'):
    sample_data[col+'_pc'] = sample_data[col]/sample_data['totalMatching']
# split into two panes assuming about 100 samples per batch
split_idx = sample_data.shape[0]//2
d1 = sample_data.iloc[:split_idx, :]
d2 = sample_data.iloc[split_idx:, :]
# plot - size optimized for 96 samples
fig, axs = plt.subplots(2,1,figsize=(12, 8))
for (i,d) in enumerate([d1, d2]):
    axs[i] = d[['final_pc',
                'failedClustering_pc',
                'failedQuality_pc',
                'failedPairProcessing_pc',
                'failedMinLen_pc',
                'failedMaxLen_pc',
                'failedNs_pc',
                'failedPossibleContamination_pc']] \
                    .plot(kind='bar', 
                          ax=axs[i],
                          legend=(True if i==0 else False),
                          stacked=True, 
                          title=('Read status breakdown per sample' if i==0 else ''))
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "filter_per_sample.pdf"), dpi=150)

# text summary 
# read counts
low_yield_cutoff = sample_data.totalMatching.mean() * PROP_READ_COUNT
low_yield_samples = sample_data.loc[sample_data.totalMatching <= low_yield_cutoff,
                                    ['totalMatching']]
msg('Samples with low number of reads matching primers. Cutoff: {}'.format(low_yield_cutoff))
msg(low_yield_samples.to_string() + '\n')
# filtering 
low_retained_samples = sample_data.loc[sample_data.final_pc <= PROP_RETAINED,
                                       ['totalMatching', 'final', 'final_pc']]
msg('Samples with at most {} reads retained after filtering'.format(PROP_RETAINED))
msg(low_retained_samples.to_string() + '\n')

# allele counts (TODO - discrete color bar)
allele_counts = pop_data.groupby(['s_Sample', 'target'], as_index=False).count()\
    .pivot(index='target', columns='s_Sample', values='h_popUID') \
    .fillna(0)
max_allele_count = allele_counts.max().max().astype(int) + 1
plot_heatmap(fill_sample_target(allele_counts),
             title="Alleles per sample per target",
             fname="allele_counts.pdf",
             cmap="coolwarm",
             center=2,
             cbar_kws=dict(ticks=range(max_allele_count)))

# major allele frequency
major_hap_freq = pop_data.groupby(['s_Sample', 'target'], as_index=False)['c_AveragedFrac'].max()\
    .pivot(index='target', columns='s_Sample', values='c_AveragedFrac')
# ignore perfectly homozygous samples
major_hap_freq = major_hap_freq.replace(1, np.nan)
# plot only if any heterozygous sites exist
if ~major_hap_freq.isna().all(axis=None):
    plot_heatmap(fill_sample_target(major_hap_freq), 
                 title="Major allele imbalance: red - over 0.5, blue - below 0.5",
                 fname="allele_imbalance.pdf",
                 cmap="coolwarm",
                 center=0.5)
else:
    print('No heterozygous calls')
    
# allele frequencies versus read counts
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.scatter(x=pop_data.c_AveragedFrac, y=np.log10(pop_data.c_ReadCnt.astype(float)), alpha=0.1)
plt.xlabel('Allele fraction')
plt.ylabel('Read count, log10')
plt.savefig(os.path.join(OUT_DIR, "allele_freq_cov.pdf"), dpi=150)

# Example:
# snakemake --use-conda --profile lsf -s blast_nt.smk/

workdir: '../../../data/phylo_ampl_dada2/coi_its2/work/'

rule all:
    input: expand('blast_xml/plate{i}.xml', i=range(1,5)) 

rule blast_nt_sanger:
    input: q='seqman_fa/{plate}.fas'
    params: d='/data/blastdb/Supported/NT/nt'
    output: 'blast_xml/{plate}.xml'
    conda: 'blast_nt.yml'
    shell: 'blastn -query {input.q} -db {params.d} -outfmt 5 > {output}'
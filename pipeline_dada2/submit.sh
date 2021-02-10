#! /bin/bash

set -eu
set -o pipefail

PARENT_DIR="/../../personal/phylo_ampl_dada2/"

# positional argument:
# directory name in tracking, used to pick up samples spreadsheet
sample_set=$1
# further arguments are passed to snakemake
smk_args=${@:2}

# directory where analysis will run
work_dir=${PARENT_DIR}${sample_set}
# path to this submission script
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# path to directory with sample set files
tracking_dir=${script_dir}/../tracking/${sample_set}

# create working directory
mkdir -p ${work_dir}
mkdir -p ${work_dir}/logs/

# copy pipeline to working directory
echo Copying pipeline to working directory ${work_dir}
cp -r ${script_dir} ${work_dir}

# update samples file in case of changes 
# otherwise pipeline starts from scratch every time due to timestamp change
echo Loading samples from ${tracking_dir}
if [ ! -f ${work_dir}/samples.csv ] || ( ! diff ${tracking_dir}/samples.csv ${work_dir}/samples.csv ); then
    cp -v ${tracking_dir}/samples.csv ${work_dir}/samples.csv
fi

# goto working directory
cd ${work_dir}

# LSF submission commands
SUB="bsub -n {cluster.cores} -q normal \
    -M {cluster.mem} \
    -R 'select[mem>{cluster.mem}] rusage[mem={cluster.mem}] span[hosts=1]' \
    -o ${work_dir}/logs/{rule}.%J.o -e ${work_dir}/logs/{rule}.%J.e"

# run pipeline
snakemake \
    -s pipeline_dada2/Snakefile \
    --cluster "${SUB}" \
    --cluster-config pipeline_dada2/cluster.json \
    --jobs 20 \
    --jn "{rulename}.{jobid}" \
    --latency-wait 30 \
    --local-cores 2 \
    --rerun-incomplete \
    --use-conda \
    --conda-prefix ../.snakemake/conda \
    ${smk_args}

# any files to copy to tracking?

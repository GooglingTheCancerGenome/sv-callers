#!/usr/bin/env bash

set -xe

MY_ENV=wf

eval "$(conda shell.bash hook)"
printenv
conda install -n base -c conda-forge mamba
conda activate base
conda list
mamba env create -n $MY_ENV -f environment.yaml
conda activate $MY_ENV
conda list
cd ./workflow
snakemake --use-conda --conda-create-envs-only --cores 1

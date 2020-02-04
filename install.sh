#!/usr/bin/env bash

set -xe

MY_ENV=wf

eval "$(conda shell.bash hook)"
printenv
conda env create -n $MY_ENV -f environment.yaml
conda activate $MY_ENV
conda list

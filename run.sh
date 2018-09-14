#!/bin/bash -xe

source ~/.profile
git lfs clone -b dev https://github.com/GooglingTheCancerGenome/sv-callers.git
cd sv-callers/snakemake
snakemake -C echo_run=0 enable_callers="['manta','delly','lumpy','gridss']" \
  mode=$1 --use-conda --latency-wait 30 --jobs \
  --cluster "xenon -vvv scheduler $2 --location local:// submit \
  --name smk.{rule} --procs-per-node=1 --start-single-process --inherit-env --max-run-time 15 --working-directory ."

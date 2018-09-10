#!/bin/bash

source ~/.profile && \
cd sv-callers/snakemake && \
git checkout dev && \
snakemake -C echo_run=0 enable_callers="['manta']" mode=$1 --latency-wait 30 --jobs \
  --cluster "xenon -vvvv scheduler $2 --location local:// submit \
  --name smk.{rule} --inherit-env --max-run-time 15 --working-directory ."

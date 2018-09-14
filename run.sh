#!/bin/bash

source ~/.profile && \
git lfs clone -b dev https://github.com/GooglingTheCancerGenome/sv-callers.git && \
cd sv-callers/snakemake && \
ls -lh data/bam/3 && \
ls -lh data/fasta && \
snakemake -C echo_run=0 enable_callers="['manta','delly','lumpy','gridss']" \
  mode=$1 --use-conda --latency-wait 30 --jobs \
  --cluster "xenon -vvv scheduler $2 --location local:// submit \
  --name smk.{rule} --inherit-env --max-run-time 15 --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log"

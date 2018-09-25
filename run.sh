#!/bin/bash -xe

source ~/.profile
git clone -b dev https://github.com/GooglingTheCancerGenome/sv-callers.git
cd sv-callers/snakemake

ECHO=$1
MODE=$2
SCH=$3
SAMPLES=$([ "$ECHO" -eq "1" ] && echo "samples_dummy.csv" || echo "samples.csv")
USE_CONDA=$([ "$ECHO" -eq "0" ] && echo "--use-conda" || echo "")

snakemake -C echo_run=$ECHO samples=$SAMPLES mode=$MODE \
  enable_callers="['manta','delly','lumpy','gridss']" $USE_CONDA \
  --latency-wait 30 --jobs \
  --cluster "xenon -vvv scheduler $SCH --location local:// submit \
  --name smk.{rule} --procs-per-node=1 --start-single-process --inherit-env \
  --max-run-time 15 --working-directory ."

if [ "$ECHO" -eq "0" ]; then
  find data -name workspace | xargs rm -fr
  find data -name *.vcf | xargs grep -cv "#"
fi

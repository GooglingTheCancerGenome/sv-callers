#!/bin/bash -xe

source ~/.profile
git clone -b dev https://github.com/GooglingTheCancerGenome/sv-callers.git
cd sv-callers/snakemake

CALLERS=(manta delly lumpy gridss)
STR_CALLERS="[$(printf "'%s'," "${CALLERS[@]}"|sed 's/,$//')]"
EXIT_CODE=0
ECHO=$1
MODE=$2
SCH=$3
SAMPLES=$([ "$ECHO" -eq "1" ] && echo "samples_dummy.csv" || echo "samples.csv")
USE_CONDA=$([ "$ECHO" -eq "0" ] && echo "--use-conda" || echo "")

echo "Selected callers: $STR_CALLERS"
snakemake -C echo_run=$ECHO samples=$SAMPLES mode=$MODE \
  enable_callers="$STR_CALLERS" $USE_CONDA \
  --configfile analysis_test.yaml \
  --latency-wait 60 --jobs \
  --cluster "xenon -vvv scheduler $SCH --location local:// submit \
  --name smk.{rule} --procs-per-node=1 --start-single-process --inherit-env \
  --max-run-time 15 --working-directory ."

if [ "$ECHO" -eq "0" ]; then
  for caller in "${CALLERS[@]}"; do
    VCF_FILE="$(find data -name $caller.vcf)"
    BOOL=$([ -e  "$VCF_FILE" ] && echo 0 || echo 1)
    MSG=$([ $BOOL -eq 0 ] && echo "Yes" || echo "No")
    EXIT_CODE=$([ $BOOL -gt $EXIT_CODE ] && echo 1 || echo 0)
    echo "$caller: VCF outfile $VCF_FILE...$MSG"
  done
fi
exit $EXIT_CODE

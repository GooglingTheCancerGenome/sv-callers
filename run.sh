#!/usr/bin/env bash

set -e

# check input arg(s)
if [ $# -lt "3" ]; then
  echo "Usage: $0 [ECHO run] [MODE {s,p}] [SCHEDULER {gridengine,slurm}]"
  exit 1
fi

ECHO="$1"
MODE="$2"
SCH="$3"
CALLERS=(manta delly lumpy gridss)
STR_CALLERS="[$(printf "'%s'," "${CALLERS[@]}"|sed 's/,$//')]"
THREADS=1
JOBS=() # array of job IDs
CONFIG="../config/analysis.yaml"
USE_CONDA=$([ "$ECHO" -eq "0" ] && echo "--use-conda" || echo "")
MY_ENV=wf

monitor () {  # monitor a job via Xenon CLI
  xenon --json scheduler "$SCH" --location local:// list --identifier "$1"
}

edit() {  # edit config file
  if [ "$ECHO" == "1" ]; then
    sed -i.org 's/echo_run:\s*0/echo_run: 1/' "$1"
  else
    sed -i.org 's/echo_run:\s*1/echo_run: 0/' "$1"
  fi

  if [ "$MODE" == "p" ]; then
    sed -i.org 's/mode:\s*s/mode: p/' "$CONFIG"
  else
    sed -i.org 's/mode:\s*p/mode: s/' "$CONFIG"
  fi
}

# activate conda env
eval "$(conda shell.bash hook)"
conda activate "$MY_ENV"
conda list

# edit CONFIG
cd ./workflow && ls -alh
edit "$CONFIG"

# run workflow
echo "Selected callers: $STR_CALLERS"
snakemake "$USE_CONDA" \
  --configfile "$CONFIG" \
  --latency-wait 60 --jobs \
  --cluster "xenon -vvv scheduler $SCH --location local:// submit \
  --name smk.{rule} --cores-per-task $THREADS --inherit-env \
  --max-run-time 15 --working-directory . \
  --stderr stderr-%j.log --stdout stdout-%j.log"

# check SV callers' output in VCF files
echo "-----------------"
echo "VCF output files:"
if [ "$ECHO" -eq "0" ]; then
  for caller in "${CALLERS[@]}"; do
    VCF_FILE=$(find data -mindepth 6 -name "$caller.vcf")
    BOOL=$([ -e  "$VCF_FILE" ] && echo 0 || echo 1)
    INFO=$([ "$BOOL" -eq "0" ] && echo "$VCF_FILE" || echo "None")
    echo " $caller: $INFO"
  done
  echo " merge:" $(find data -name all.vcf)
fi

# show logs
echo "----------"
echo "Log files:"
for f in *.log; do
  echo -e "\n### $f ###\n"
  cat "$f"
  JOB_ID=$(echo "$f" | grep err | sed 's/\w*-\([0-9]*\).\w*/\1/')
  JOBS+=($JOB_ID)
done

# collect job accounting info
for j in ${JOBS[@]}; do
  monitor "$j" > "jobacct-$j.json"
done
cat ./*.json

# exit with non-zero if there are failed jobs
[[ $(jq ".statuses | .[] | select(.done==true and .exitCode!=0)" ./*.json) ]] \
  && exit 1 || exit 0

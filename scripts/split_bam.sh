#!/bin/bash -e
#
# This script splits a whole-genome BAM file into per-chromosome BAM files
# on SGE cluster. Note that only proper read pairs are selected per chromosome.
#

BAM_FILE=$1
THREADS=$2

if [ ! -f ${BAM_FILE} ]; then
  echo "${BAM_FILE} does not exist."
  exit 1
fi

if [ "${THREADS}" -lt 1 ]; then
  echo "Set the number of threads to a positive integer >= 1."
  exit 1
fi

BAM_PREFIX=$(basename ${BAM_FILE} .bam)
RTIME=7200 # in sec
VMEM=100M

# fetch a list of chromosomes from the BAM header
STR=$(echo $(samtools view -H ${BAM_FILE}|grep ^@SQ|cut -f 2|cut -f 2 -d :))
IFS=' ' read -ra CHR <<< "${STR}"

if [ "${#CHR[@]}" -eq 0 ]; then
  echo "The list of chromosomes is empty."
  exit 1
fi

for c in "${CHR[@]}"; do
  CHR_BAM_FILE=${BAM_PREFIX}-chr${c}.bam
  echo "
    samtools view -bh -@ ${THREADS} -f 2 -o ${CHR_BAM_FILE} ${BAM_FILE} ${c} \
    && samtools index -@ ${THREADS} ${CHR_BAM_FILE} ${c}.bai"| \
    qsub -N samtools-${THREADS}-${CHR_BAM_FILE} -cwd -V \
    -l h_rt=$((${RTIME}/${THREADS})) -l h_vmem=${VMEM} -pe threaded ${THREADS}
done

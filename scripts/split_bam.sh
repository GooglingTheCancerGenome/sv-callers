#!/bin/bash
#
# This script splits a whole-genome BAM file into per-chromosome BAM files
# on SGE cluster. Note that only proper read pairs are selected per chromosome.
#

BAM_FILE=$1
THREADS=$2

if [ ! -f $BAM_FILE ]; then
    echo "$BAM_FILE does not exist."
    exit 1
fi

if ! [ "$THREADS" -ge 1 ]; then
   echo "Set the number of threads to a positive integer >= 1."
   exit 1
fi

# fetch list of chromosomes from the BAM header
STR=$(echo $(samtools view -H $BAM_FILE | grep ^@SQ|cut -f 2 | cut -f 2 -d :))
IFS=' ' read -ra CHR <<< "$STR"

for c in "${CHR[@]}"
do
    echo "
    samtools view -bh -@ $THREADS -f 2 -o ${c}.bam VCAP_dedup.realigned.bam $c \
    && samtools index -@ $THREADS ${c}.bam ${c}.bai
    " | qsub -N samtools-${THREADS}-chr${c} -cwd -V -l h_rt=02:00:00 -l h_vmem=100M #-pe threaded $THREADS
done

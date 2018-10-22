#!/bin/bash

# This bash script generates the SURVIVOR merge VCF file with the overlap among sv callers within each run
# SURVIVOR version v1.0.5 was used

for SAMPLE in COLO829 NA12878; do
	echo "Sample: "$SAMPLE
	for RUN in `seq 1 10`; do
		echo "Run: "$RUN
		INPUT_FILES="/hpc/cog_bioinf/ridder/users/akuzniar/"$SAMPLE"/bam/"$RUN"/"*"/"*"/"[a-z]*".vcf"
		OUTPUT_DIR="/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Papers/SV_callers/"$SAMPLE"/all/"
		VCF_FILES=$OUTPUT_DIR"/"$SAMPLE"_run"$RUN"_VCFs.txt"
		OUTPUT_FILE=$OUTPUT_DIR"/"$SAMPLE"_run"$RUN"_SURVIVOR.vcf"
		ls $INPUT_FILES | grep -v 'unfiltered' > $VCF_FILES
		SURVIVOR merge $VCF_FILES 100 1 0 0 0 0 $OUTPUT_FILE
	done
done

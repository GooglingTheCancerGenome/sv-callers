#!/bin/bash

# This bash script generates the SURVIVOR merge VCF file with the overlap among runs for each sv caller
# SURVIVOR version v1.0.5 was used

for SAMPLE in COLO829 NA12878; do
	echo "Sample: "$SAMPLE
	for CALLER in delly lumpy manta gridss; do
		echo "Caller: "$CALLER
		INPUT_FILES="/hpc/cog_bioinf/ridder/users/akuzniar/"$SAMPLE"/bam/"*"/"*"/"$CALLER*"/"$CALLER".vcf"
		OUTPUT_DIR="/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Papers/SV_callers/"$SAMPLE"/"$CALLER"/"
		VCF_FILES=$OUTPUT_DIR"/"$SAMPLE"_"$CALLER"_VCFs.txt"
		OUTPUT_FILE=$OUTPUT_DIR"/"$SAMPLE"_"$CALLER"_SURVIVOR.vcf"
		ls $INPUT_FILES > $VCF_FILES
		SURVIVOR merge $VCF_FILES 100 1 0 0 0 0 $OUTPUT_FILE
	done
done

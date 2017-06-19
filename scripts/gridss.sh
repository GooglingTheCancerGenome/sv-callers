#!/bin/bash

#       -Dsamjdk.use_async_io_read_samtools=true \
#       -Dsamjdk.use_async_io_write_samtools=true \
#       -Dsamjdk.use_async_io_write_tribble=true \
#       -Dsamjdk.compression_level=1 \

gridss	-Xmx16g \
	gridss.CallVariants \
	WORKER_THREADS=8 \
	TMP_DIR=. \
	WORKING_DIR=. \
	REFERENCE_SEQUENCE="Homo_sapiens.GRCh37.GATK.illumina.fasta" \
	INPUT="chr21_pairs.bam" \
	OUTPUT="gridss_chr21.vcf" \
	ASSEMBLY="gridss_assembly.bam"

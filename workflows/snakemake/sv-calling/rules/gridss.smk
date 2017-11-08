rule gridss:
    input:
        REF_FASTA,
        REF_FAI,
        TUMOR_BAM,
        TUMOR_BAI,
        #TODO: Add BWA 0.6.x index files (.amb, .ann, .bwt, .pac, .sa)
    output:
        dir = "gridss_out",
        vcf = "gridss_out/%s.vcf" % TUMOR
    conda:
        "envs/sv_callers.yaml"
    threads: 8
    shell:
        """
        gridss -Xmx31g gridss.CallVariants WORKER_THREADS={threads} \
            REFERENCE_SEQUENCE={REF_FASTA} \
            INPUT={TUMOR_BAM} \
            OUTPUT={output.vcf} \
            ASSEMBLY={TUMOR}_assembly.bam \
            WORKING_DIR={output.dir} \
            TMP_DIR={output.dir}
        """

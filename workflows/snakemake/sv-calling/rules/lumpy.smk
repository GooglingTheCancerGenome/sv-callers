rule lumpy:
    input:
        REF_FASTA,
        REF_FAI,
        TUMOR_BAM,
        TUMOR_BAI,
    output:
        #dir = "lumpy_out",
        vcf = "lumpy_out/%s.vcf" % TUMOR
    threads: 1
    conda:
        "envs/sv_callers.yaml"
    shell:
        """
        lumpyexpress -B {TUMOR_BAM} -o {output.vcf}
        """

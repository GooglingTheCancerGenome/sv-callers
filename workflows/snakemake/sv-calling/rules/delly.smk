rule delly:
    input:
        REF_FASTA,
        REF_FAI,
        TUMOR_BAM,
        TUMOR_BAI,
        excl = TUMOR + ".excl"
    output:
        #dir = "delly_out",
        bcf = "delly_out/%s-%s.bcf" % (TUMOR, SV_TYPE)
    conda:
        "envs/sv_callers.yaml"
    threads: 2
    shell:
        """
        delly call -t {SV_TYPE} \
            -g {REF_FASTA} \
            -x {input.excl} \
            -o {output.bcf} {TUMOR_BAM}
        """

rule bcf_to_vcf:
    input:
        "delly_out/%s-%s.bcf" % (TUMOR, SV_TYPE)
    output:
        "delly_out/%s-%s.vcf" % (TUMOR, SV_TYPE)
    conda:
      "envs/sv_callers.yaml"
    shell:
      "bcftools view {input} > {output}"

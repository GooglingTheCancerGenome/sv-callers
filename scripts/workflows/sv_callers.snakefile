REF = "Homo_sapiens.GRCh37.GATK.illumina"
REF_FASTA = REF + ".fasta"
REF_FAI = REF_FASTA + ".fai"
SV_TYPE = "BND"
TUMOR = "chr22"
TUMOR_BAM = TUMOR + ".bam"
TUMOR_BAI = TUMOR + ".bai"

rule manta:
    input:
        REF_FASTA,
        REF_FAI,
        TUMOR_BAM,
        TUMOR_BAI,
    output:
        dir = "manta_out"
    conda:
        "envs/sv_callers.yaml"
    threads: 16
    shell:
        """
        configManta.py --tumorBam {TUMOR_BAM} \
            --reference {REF_FASTA} \
            --runDir {output.dir}
        cd {output.dir} && ./runWorkflow.py --quiet -m local -j {threads}
        """

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
    threads: 1
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


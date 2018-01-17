rule gridss:
    input:
        fasta = get_fasta(),
        fai = get_faidx(),  # bwa index files also required
        tumor_bam = "{sampledir}/{tumor}" + get_filext("bam"),
        tumor_bai = "{sampledir}/{tumor}" + get_filext("bam_idx"),
        normal_bam = "{sampledir}/{normal}" + get_filext("bam"),
        normal_bai = "{sampledir}/{normal}" + get_filext("bam_idx")
    params:
        outdir = os.path.join("{sampledir}", get_outdir("gridss"))
    output:
        log = os.path.join("{sampledir}", get_outdir("gridss"),
                           "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("gridss")
    resources:
        mem_mb = get_maxmem("gridss")
    shell:
        """
        if [ "{config[echo_run]}" = "1" ]; then
            echo "{input}" > "{output}"
        else
            # clean-up prior to SV calling
            rm -f "{input.fasta}.dict" && \
            gridss -Xmx31g gridss.CallVariants WORKER_THREADS={threads} \
                REFERENCE_SEQUENCE="{input.fasta}" \
                INPUT="{input.normal_bam}" \
                INPUT="{input.tumor_bam}" \
                OUTPUT="{params}/gridss.vcf" \
                ASSEMBLY="{params}/assembly.bam" \
                WORKING_DIR="{params}" \
                TMP_DIR="{params}" 2>&1
        fi
        """

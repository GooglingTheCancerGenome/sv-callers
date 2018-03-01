rule gridss:
    input:
        fasta = get_fasta(),
        fai = get_faidx(),  # bwa index files also required
        tumor_bam = "{path}/{tumor}" + get_filext("bam"),
        tumor_bai = "{path}/{tumor}" + get_filext("bam_idx"),
        normal_bam = "{path}/{normal}" + get_filext("bam"),
        normal_bai = "{path}/{normal}" + get_filext("bam_idx")
    output:
        log = os.path.join("{path}", "{tumor}--{normal}", get_outdir("gridss"),
                           "{rule}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("gridss")
    resources:
        mem_mb = get_memory("gridss"),
        tmp_mb = get_tmpspace("gridss")
    shell:
        """
        # if 'tmpspace' set to >0MB use TMPDIR otherwise use OUTDIR
        OUTDIR="$(dirname "{output}")"
        TMP=$([ "{resources.tmp_mb}" = "0" ] &&
            echo "${{OUTDIR}}" ||
            echo "${{TMPDIR}}")
        if [ "{config[echo_run]}" = "1" ]; then
            echo "{input}" "${{TMP}}" > "{output}"
        else
            # clean-up prior to SV calling
            rm -f "{input.fasta}.dict" &&
            gridss -Xmx31g gridss.CallVariants \
                WORKER_THREADS={threads} \
                REFERENCE_SEQUENCE="{input.fasta}" \
                INPUT="{input.normal_bam}" \
                INPUT="{input.tumor_bam}" \
                OUTPUT="${{OUTDIR}}/gridss.vcf" \
                ASSEMBLY="${{OUTDIR}}/gridss_assembly.bam" \
                WORKING_DIR="${{OUTDIR}}" \
                TMP_DIR="${{TMP}}/gridss.${{RANDOM}}" 2>&1
            date "+%Y-%m-%d %H:%M:%S" > "{output}"
        fi
        """

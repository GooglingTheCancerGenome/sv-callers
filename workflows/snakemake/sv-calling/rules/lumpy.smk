rule lumpy:
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = "{sampledir}/{tumor}" + get_filext("bam"),
        tumor_bai = "{sampledir}/{tumor}" + get_filext("bam_idx"),
        normal_bam = "{sampledir}/{normal}" + get_filext("bam"),
        normal_bai = "{sampledir}/{normal}" + get_filext("bam_idx")
    params:
        outdir = os.path.join("{sampledir}", get_outdir("lumpy"))
    output:
        log = os.path.join("{sampledir}", get_outdir("lumpy"),
                           "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("lumpy")
    resources:
        mem_mb = get_memory("lumpy"),
        tmp_mb = get_tmpspace("lumpy")
    shell:
        """
        # TMPDIR used only if 'tmp' is >0 MB in config
        TMP=$([ "{resources.tmp_mb}" = "0" ] && echo "{params}" ||
            echo "${{TMPDIR}}")
        if [ "{config[echo_run]}" = "1" ]; then
            echo "{input}" "${{TMP}}"> "{output}"
        else
            lumpyexpress \
                -B "{input.tumor_bam}","{input.normal_bam}" \
                -o "{params}/lumpy.vcf" \
                -m 4 `# min. sample weight` \
                -r 0 `# trim threshold` \
                -T "${{TMP}}" 2>&1
            date "+%Y-%m-%d %H:%M:%S" > "{output}"
        fi
        """

rule lumpy_s:  # somatic mode
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = "{path}/{tumor}" + get_filext("bam"),
        tumor_bai = "{path}/{tumor}" + get_filext("bam_idx"),
        normal_bam = "{path}/{normal}" + get_filext("bam"),
        normal_bai = "{path}/{normal}" + get_filext("bam_idx")
    params:
        excl_opt = '-x "%s"' % get_bed("lumpy") if get_bed("lumpy") else ""
    output:
        os.path.join("{path}/{tumor}--{normal}", get_outdir("lumpy"),
                     "lumpy" + get_filext("vcf"))
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("lumpy")
    resources:
        mem_mb = get_memory("lumpy"),
        tmp_mb = get_tmpspace("lumpy")
    shell:
        """
        set -x

        # if 'tmpspace' set to >0MB use TMPDIR otherwise use OUTDIR
        OUTDIR="$(dirname "{output}")"
        TMP=$([ "{resources.tmp_mb}" -eq "0" ] &&
            echo "${{OUTDIR}}" || echo "${{TMPDIR}}")

        # run dummy or real job
        if [ "{config[echo_run]}" -eq "1" ]; then
            echo "{input}" "${{TMP}}" > "{output}"
        else
            lumpyexpress \
                -B "{input.tumor_bam}","{input.normal_bam}" \
                {params.excl_opt} \
                -o "{output}" \
                -m 4 `# min. sample weight` \
                -r 0 `# trim threshold` \
                -k `# keep tmp files` \
                -T "${{TMP}}/lumpy.${{RANDOM}}"
        fi
        """

rule lumpy_g:  # germline mode
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        bam = "{path}/{sample}" + get_filext("bam"),
        bai = "{path}/{sample}" + get_filext("bam_idx")
    params:
        excl_opt = '-x "%s"' % get_bed("lumpy") if get_bed("lumpy") else ""
    output:
        os.path.join("{path}/{sample}", get_outdir("lumpy"), "lumpy" +
                     get_filext("vcf"))
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("lumpy")
    resources:
        mem_mb = get_memory("lumpy"),
        tmp_mb = get_tmpspace("lumpy")
    shell:
        """
        set -x

        # if 'tmpspace' set to >0MB use TMPDIR otherwise use OUTDIR
        OUTDIR="$(dirname "{output}")"
        TMP=$([ "{resources.tmp_mb}" -eq "0" ] &&
            echo "${{OUTDIR}}" || echo "${{TMPDIR}}")

        # run dummy or real job
        if [ "{config[echo_run]}" -eq "1" ]; then
            echo "{input}" "${{TMP}}" > "{output}"
        else
            lumpyexpress \
                -B "{input.bam}" \
                {params.excl_opt} \
                -o "{output}" \
                -m 4 `# min. sample weight` \
                -r 0 `# trim threshold` \
                -k `# keep tmp files` \
                -T "${{TMP}}/lumpy.${{RANDOM}}"
        fi
        """

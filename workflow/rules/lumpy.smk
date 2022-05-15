rule lumpy_p:  # paired-samples analysis
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = get_bam("{path}/{tumor}"),
        tumor_bai = get_bai("{path}/{tumor}"),
        normal_bam = get_bam("{path}/{normal}"),
        normal_bai = get_bai("{path}/{normal}")
    params:
        excl_opt = '-x "%s"' % get_bed() if exclude_regions() else ""
    output:
        os.path.join("{path}/{tumor}--{normal}", get_outdir("lumpy"),
                     "lumpy" + config.file_exts.vcf)
    conda:
        "../envs/caller.yaml"
    threads:
        config.callers.lumpy.threads
    resources:
        mem_mb = config.callers.lumpy.memory,
        tmp_mb = config.callers.lumpy.tmpspace
    shell:
        """
        set -x

        # if 'tmpspace' set to >0MB use TMPDIR otherwise use OUTDIR
        OUTDIR="$(dirname "{output}")"
        PREFIX="$(basename "{output}" .vcf)"
        OUTFILE="${{OUTDIR}}/${{PREFIX}}.unfiltered.vcf"
        TMP=$([ "{resources.tmp_mb}" -eq "0" ] &&
            echo "${{OUTDIR}}" || echo "${{TMPDIR}}")

        # run dummy or real job
        if [ "{config.echo_run}" -eq "1" ]; then
            echo "{input}" "${{TMP}}" > "{output}"
        else
            lumpyexpress \
                -B "{input.tumor_bam}","{input.normal_bam}" \
                {params.excl_opt} \
                -o "${{OUTFILE}}" \
                -m 4 `# min. sample weight` \
                -r 0 `# trim threshold` \
                -k `# keep tmp files` \
                -T "${{TMP}}/lumpy.${{RANDOM}}"
            # somatic + SV quality filtering
            #   'normal' sample assumes index 1
            bcftools filter \
                -O v `# uncompressed VCF format` \
                -o "{output}" \
                -i "FORMAT/SU[1] == 0 && FILTER == '.'" \
                "${{OUTFILE}}"
        fi
        """

rule lumpy_s:  # single-sample analysis
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        bam = get_bam("{path}/{sample}"),
        bai = get_bai("{path}/{sample}")
    params:
        excl_opt = '-x "%s"' % get_bed() if exclude_regions() else ""
    output:
        os.path.join("{path}/{sample}", get_outdir("lumpy"), "lumpy" +
                     config.file_exts.vcf)
    conda:
        "../envs/caller.yaml"
    threads:
        config.callers.lumpy.threads
    resources:
        mem_mb = config.callers.lumpy.memory,
        tmp_mb = config.callers.lumpy.tmpspace
    shell:
        """
        set -x

        # if 'tmpspace' set to >0MB use TMPDIR otherwise use OUTDIR
        OUTDIR="$(dirname "{output}")"
        PREFIX="$(basename "{output}" .vcf)"
        OUTFILE="${{OUTDIR}}/${{PREFIX}}.unfiltered.vcf"
        TMP=$([ "{resources.tmp_mb}" -eq "0" ] &&
            echo "${{OUTDIR}}" || echo "${{TMPDIR}}")

        # run dummy or real job
        if [ "{config.echo_run}" -eq "1" ]; then
            echo "{input}" "${{TMP}}" > "{output}"
        else
            lumpyexpress \
                -B "{input.bam}" \
                {params.excl_opt} \
                -o "${{OUTFILE}}" \
                -m 4 `# min. sample weight` \
                -r 0 `# trim threshold` \
                -k `# keep tmp files` \
                -T "${{TMP}}/lumpy.${{RANDOM}}"
            # SV quality filtering
            bcftools filter \
                -O v `# uncompressed VCF format` \
                -o "{output}" \
                -i "FILTER == '.'" \
                "${{OUTFILE}}"
        fi
        """

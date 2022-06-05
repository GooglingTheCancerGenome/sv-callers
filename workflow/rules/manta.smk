rule manta_p:  # paired-samples analysis
    input:
        fasta=get_fasta(),
        fai=get_faidx()[0],
        tumor_bam=get_bam("{path}/{tumor}"),
        tumor_bai=get_bai("{path}/{tumor}"),
        normal_bam=get_bam("{path}/{normal}"),
        normal_bai=get_bai("{path}/{normal}"),
    params:
        #    excl_opt = '-x "%s"' % get_bed() if exclude_regions() else ""
        outfile="results/variants/somaticSV.vcf.gz",
    output:
        os.path.join(
            "{path}/{tumor}--{normal}",
            get_outdir("manta"),
            "manta{}".format(config.file_exts.vcf),
        ),
    conda:
        "../envs/caller.yaml"
    threads: config.callers.manta.threads
    resources:
        mem_mb=config.callers.manta.memory,
        tmp_mb=config.callers.manta.tmpspace,
    shell:
        """
        set -xe

        OUTDIR="$(dirname "{output}")"
        OUTFILE="$(basename "{output}")"

        # run dummy or real job
        if [ "{config.echo_run}" -eq "1" ]; then
            echo "{input}" > "{output}"
        else
            configManta.py \
                --runDir "${{OUTDIR}}" \
                --reference "{input.fasta}" \
                --tumorBam "{input.tumor_bam}" \
                --normalBam "{input.normal_bam}"
            cd "${{OUTDIR}}"
            ./runWorkflow.py \
                --quiet \
                -m local \
                -j {threads}
            # SV quality filtering
            bcftools filter \
                -O v `# uncompressed VCF format` \
                -o "${{OUTFILE}}" \
                -i "FILTER == 'PASS'" \
                "{params.outfile}"
        fi
        """


rule manta_s:  # single-sample analysis: germline or tumor-only
    input:
        fasta=get_fasta(),
        fai=get_faidx()[0],
        bam=get_bam("{path}/{sample}"),
        bai=get_bai("{path}/{sample}"),
    params:
        # excl_opt = '-x "%s"' % get_bed() if exclude_regions() else ""
        bam_opt="--tumorBam" if is_tumor_only() else "--bam",
        outfile="results/variants/tumorSV.vcf.gz"
        if is_tumor_only()
        else "results/variants/diploidSV.vcf.gz",
    output:
        os.path.join(
            "{path}/{sample}",
            get_outdir("manta"),
            "manta{}".format(config.file_exts.vcf),
        ),
    conda:
        "../envs/caller.yaml"
    threads: config.callers.manta.threads
    resources:
        mem_mb=config.callers.manta.memory,
        tmp_mb=config.callers.manta.tmpspace,
    shell:
        """
        set -xe

        OUTDIR="$(dirname "{output}")"
        OUTFILE="$(basename "{output}")"

        # run dummy or real job
        if [ "{config.echo_run}" -eq "1" ]; then
            echo "{input}" > "{output}"
        else
            configManta.py \
                --runDir "${{OUTDIR}}" \
                --reference "{input.fasta}" \
                {params.bam_opt} "{input.bam}"
            cd "${{OUTDIR}}"
            ./runWorkflow.py \
                --quiet \
                -m local \
                -j {threads}
            # SV quality filtering
            bcftools filter \
                -O v `# uncompressed VCF format` \
                -o "${{OUTFILE}}" \
                -i "FILTER == 'PASS'" \
                "{params.outfile}"
        fi
        """

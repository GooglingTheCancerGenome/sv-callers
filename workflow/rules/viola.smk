rule viola:  # used by both modes
    input:
        os.path.join(
            "{path}",
            "{tumor}--{normal}",
            "{outdir}",
            "{}{}".format("{prefix}", config.file_exts.vcf),
        )
        if config.mode.PAIRED_SAMPLE is True
        else os.path.join(
            "{path}",
            "{sample}",
            "{outdir}",
            "{}{}".format("{prefix}", config.file_exts.vcf),
        ),
    output:
        os.path.join(
            "{path}",
            "{tumor}--{normal}",
            "{outdir}",
            "viola",
            "{}{}".format("{prefix}", config.file_exts.vcf),
        )
        if config.mode.PAIRED_SAMPLE is True
        else os.path.join(
            "{path}",
            "{sample}",
            "{outdir}",
            "viola",
            "{}{}".format("{prefix}", config.file_exts.vcf),
        ),
    conda:
        "../envs/postproc.yaml"
    threads: config.postproc.survivor.threads
    resources:
        mem_mb=config.postproc.survivor.memory,
        tmp_mb=config.postproc.survivor.tmpspace,
    script:
        "../scripts/viola_vcf.py"

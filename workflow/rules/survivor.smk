rule survivor_filter:  # used by both modes
    input:
        os.path.join(
            "{path}",
            "{tumor}--{normal}",
            "{outdir}",
            "viola",
            "{}{}".format("{prefix}", config.file_exts.vcf),
        )
        if config.mode is config.mode.PAIRED_SAMPLE
        else os.path.join(
            "{path}",
            "{sample}",
            "{outdir}",
            "viola",
            "{}{}".format("{prefix}", config.file_exts.vcf),
        ),
    output:
        os.path.join(
            "{path}",
            "{tumor}--{normal}",
            "{outdir}",
            get_outdir("survivor"),
            "{}{}".format("{prefix}", config.file_exts.vcf),
        )
        if config.mode is config.mode.PAIRED_SAMPLE
        else os.path.join(
            "{path}",
            "{sample}",
            "{outdir}",
            get_outdir("survivor"),
            "{}{}".format("{prefix}", config.file_exts.vcf),
        ),
    params:
        excl=exclude_regions(),
        args=survivor_args("filter"),
    conda:
        "../envs/postproc.yaml"
    threads: config.postproc.survivor.threads
    resources:
        mem_mb=config.postproc.survivor.memory,
        tmp_mb=config.postproc.survivor.tmpspace,
    shell:
        """
        set -xe

        # run dummy or real job
        if [ "{config.echo_run}" -eq "1" ]; then
            cat "{input}" > "{output}"
        else
            if [ "{params.excl}" -eq "1" ]; then
                SURVIVOR filter "{input}" {params.args} "{output}"
            else
                ln -sr "{input}" "{output}"
            fi
        fi
        """


rule survivor_merge:  # used by both modes
    input:
        [
            os.path.join(
                "{path}",
                "{tumor}--{normal}",
                get_outdir(c),
                get_outdir("survivor"),
                "{}{}".format(c, config.file_exts.vcf),
            )
            for c in config.enable_callers
        ]
        if config.mode is config.mode.PAIRED_SAMPLE
        else [
            os.path.join(
                "{path}",
                "{sample}",
                get_outdir(c),
                get_outdir("survivor"),
                "{}{}".format(c, config.file_exts.vcf),
            )
            for c in config.enable_callers
        ],
    params:
        args=survivor_args("merge")[1:-1],
    output:
        [
            os.path.join("{path}", "{tumor}--{normal}", survivor_args("merge")[0]),
            os.path.join("{path}", "{tumor}--{normal}", survivor_args("merge")[-1]),
        ]
        if config.mode is config.mode.PAIRED_SAMPLE
        else [
            os.path.join("{path}", "{sample}", survivor_args("merge")[0]),
            os.path.join("{path}", "{sample}", survivor_args("merge")[-1]),
        ],
    conda:
        "../envs/postproc.yaml"
    threads: config.postproc.survivor.threads
    resources:
        mem_mb=config.postproc.survivor.memory,
        tmp_mb=config.postproc.survivor.tmpspace,
    shell:
        """
        set -xe

        # create a list of VCF files
        for f in $(echo "{input}")
        do
            echo "$f" >> "{output[0]}"
        done

        # run dummy or real job
        if [ "{config.echo_run}" -eq "1" ]; then
            cat "{output[0]}" > "{output[1]}"
        else
            SURVIVOR merge "{output[0]}" {params.args} "{output[1]}"
        fi
        """

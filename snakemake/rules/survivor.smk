rule survivor_filter:  # used by both modes
    input:
        os.path.join("{path}", "{tumor}--{normal}", "{outdir}", "{prefix}" + get_filext("vcf"))
        if config["mode"].startswith("p") is True else
        os.path.join("{path}", "{sample}", "{outdir}", "{prefix}" + get_filext("vcf"))
    output:
        os.path.join("{path}", "{tumor}--{normal}", "{outdir}", "survivor", "{prefix}" + get_filext("vcf"))
        if config["mode"].startswith("p") is True else
        os.path.join("{path}", "{sample}", "{outdir}", "survivor", "{prefix}" + get_filext("vcf"))
    params:
        bed = '"%s"' % get_bed() if exclude_regions() else ""
    conda:
        "../environment.yaml"
    threads: 1
    resources:
        mem_mb = 1024,
        tmp_mb = 0
    shell:
        """
        set -x

        # run dummy or real job
        if [ "{config[echo_run]}" -eq "1" ]; then
            cat "{input}" > "{output}"
        else
            if [ "{params.bed}" == "" ]; then
                ln -s "{input}" "{output}"
            else
                SURVIVOR filter "{input}" "{params.bed}" -1 -1 0 -1 "{output}"
            fi
        fi
        """

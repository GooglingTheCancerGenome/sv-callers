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
                ln -sr "{input}" "{output}"
            else
                SURVIVOR filter "{input}" "{params.bed}" -1 -1 0 -1 "{output}"
            fi
        fi
        """

rule survivor_merge:  # used by both modes
    input:
        [os.path.join("{path}", "{tumor}--{normal}", get_outdir(c), "survivor", c + get_filext("vcf"))
         for c in get_callers()]
        if config["mode"].startswith("p") is True else
        [os.path.join("{path}", "{sample}", get_outdir(c), "survivor", c + get_filext("vcf"))
         for c in get_callers()]
    output:
        os.path.join("{path}", "{tumor}--{normal}",  "all" + get_filext("vcf"))
        if config["mode"].startswith("p") is True else
        os.path.join("{path}", "{sample}", "all" + get_filext("vcf"))
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
            echo "{input}" | tr ' ' '\n' > all.vcf &&
            SURVIVOR merge all.txt 100 1 0 0 0 0 "{output}"
        fi
        """

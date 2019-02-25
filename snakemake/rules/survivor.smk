rule bed_filter:  # used by both modes
    input:
        [os.path.join("{path}/{tumor}--{normal}", get_outdir(c), c +
            get_filext("vcf")) for c in get_callers()]
        if config["mode"].startswith("p") is True else
        [os.path.join("{path}/{sample}", get_outdir(c), c +
            get_filext("vcf")) for c in get_callers()]
    output:
        [os.path.join("{path}/{tumor}--{normal}", get_outdir(c), c +
            get_filext("bed") + get_filext("vcf")) for c in get_callers()]
        if config["mode"].startswith("p") is True else
        [os.path.join("{path}/{sample}", get_outdir(c), c + get_filext("bed") +
            get_filext("vcf")) for c in get_callers()]
    params:
        bed = config["exclusion_list"]
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
            SURVIVOR filter {input} {params.bed} -1 -1 0 -1 {output}
       fi
       """

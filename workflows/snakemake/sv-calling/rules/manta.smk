rule manta:
    input:
        REF_FASTA,
        REF_FAI,
        TUMOR_BAM,
        TUMOR_BAI,
    output:
        dir = "manta_out"
    conda:
        "envs/sv_callers.yaml"
    threads: 16
    shell:
        """
        configManta.py --tumorBam {TUMOR_BAM} \
            --reference {REF_FASTA} \
            --runDir {output.dir}
        cd {output.dir} && ./runWorkflow.py --quiet -m local -j {threads}
        """

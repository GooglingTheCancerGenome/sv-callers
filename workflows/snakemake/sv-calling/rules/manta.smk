rule manta:
    input:
        fasta=config["ref_genome"] + config["file_exts"]["fasta"],
        fai=config["ref_genome"] + config["file_exts"]["fai"],
        tumor_bam="{tumor}" + config["file_exts"]["bam"],
        tumor_bai="{tumor}" + config["file_exts"]["bai"],
        normal_bam="{normal}" + config["file_exts"]["bam"],
        normal_bai="{normal}" + config["file_exts"]["bai"]
    output:
        log="{tumor}-{normal}.log"
    params:
        outdir=config["outdirs"]["manta"]
    conda:
        "../environment.yaml"
    threads: 8
    # shell:
    #     """
    #     echo {params.outdir} {input} > {output}
    #     """
    shell:
        """
        configManta.py --runDir {params.outdir} \
            --reference {input.fasta} \
            --tumorBam {input.tumor_bam} \
            --normalBam {input.normal_bam}
        cd {params.outdir} && ./runWorkflow.py --quiet -m local -j {threads}
        """

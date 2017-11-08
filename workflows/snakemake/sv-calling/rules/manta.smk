rule manta:
    input:
        fasta=config["genome"] + config["file_exts"]["fasta"],
        fai=config["genome"] + config["file_exts"]["fai"],
        tumor_bam="{path}/{tumor}" + config["file_exts"]["bam"],
        tumor_bai="{path}/{tumor}" + config["file_exts"]["bai"],
        normal_bam="{path}/{normal}" + config["file_exts"]["bam"],
        normal_bai="{path}/{normal}" + config["file_exts"]["bai"]
    params:
        config["sv_dirs"]["manta"]
    output:
        "{path}/{params}/{tumor}-{normal}.log"
    conda:
        "../environment.yaml"
    threads: 8
    # shell:
    #     """
    #     echo {input} > {output}
    #     """
    shell:
        """
        configManta.py --runDir "{path}/{params}" \
            --reference {input.fasta} \
            --tumorBam {input.tumor_bam} \
            --normalBam {input.normal_bam}
        cd {params} && ./runWorkflow.py --quiet -m local -j {threads}
        """

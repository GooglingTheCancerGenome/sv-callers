rule manta:
    input:
        fasta=config["ref_genome"] + config["file_exts"]["fasta"],
        fai=config["ref_genome"] + config["file_exts"]["fai"],
        tumor_bam="{base_dir}/{tumor}" + config["file_exts"]["bam"],
        tumor_bai="{base_dir}/{tumor}" + config["file_exts"]["bai"],
        normal_bam="{base_dir}/{normal}" + config["file_exts"]["bam"],
        normal_bai="{base_dir}/{normal}" + config["file_exts"]["bai"]
    output:
        os.path.join("{base_dir}", get_outdir("manta"), \
            "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        config["sv_callers"]["manta"]["threads"]
    shell:
        """
        echo {input} > {output}
        """
#    shell:
#        """
#        configManta.py --runDir "{wildcards.base_dir}" \
#            --reference "{input.fasta}" \
#            --tumorBam "{input.tumor_bam}" \
#            --normalBam "{input.normal_bam}"
#        cd "{wildcards.base_dir}" && ./runWorkflow.py --quiet -m local -j {threads}
#        date "+%Y-%m-%d %H:%M:%S" > "{output}"
#        """

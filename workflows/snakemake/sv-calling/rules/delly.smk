rule delly:
    input:
        fasta=config["ref_genome"] + config["file_exts"]["fasta"],
        fai=config["ref_genome"] + config["file_exts"]["fai"],
        tumor_bam="{base_dir}/{tumor}" + config["file_exts"]["bam"],
        tumor_bai="{base_dir}/{tumor}" + config["file_exts"]["bai"],
        normal_bam="{base_dir}/{normal}" + config["file_exts"]["bam"],
        normal_bai="{base_dir}/{normal}" + config["file_exts"]["bai"]
    output:
        os.path.join("{base_dir}", get_outdir("delly"), \
            "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        config["sv_callers"]["delly"]["threads"]
    resources:
        mem_mb=config["sv_callers"]["delly"]["memory"]
    shell:
        """
        echo {input} > {output}
        """
#    shell:
#        """
#        delly call -t BND \
#            -g {input.fasta} \
#            -o {wildcards.base_dir}/{params}/BND.bcf \
#            {input.tumor_bam} {input.normal_bam}
#        date "+%Y-%m-%d %H:%M:%S" > {output}
#        """

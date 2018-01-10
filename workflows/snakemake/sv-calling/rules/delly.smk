rule delly:
    input:
        fasta=config["ref_genome"] + config["file_exts"]["fasta"],
        fai=config["ref_genome"] + config["file_exts"]["fai"],
        tumor_bam="{sampledir}/{tumor}" + config["file_exts"]["bam"],
        tumor_bai="{sampledir}/{tumor}" + config["file_exts"]["bai"],
        normal_bam="{sampledir}/{normal}" + config["file_exts"]["bam"],
        normal_bai="{sampledir}/{normal}" + config["file_exts"]["bai"]
    output:
        os.path.join("{sampledir}", get_outdir("delly"), \
            "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("delly")
    resources:
        mem_mb=get_maxmem("delly")
    shell:
        """
        echo {input} > {output}
        """
#    shell:
#        """
#        delly call -t BND \
#            -g {input.fasta} \
#            -o {wildcards.sampledir}/{params}/BND.bcf \
#            {input.tumor_bam} {input.normal_bam}
#        date "+%Y-%m-%d %H:%M:%S" > {output}
#        """

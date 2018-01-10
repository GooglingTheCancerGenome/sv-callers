rule lumpy:
    input:
        fasta=config["ref_genome"] + config["file_exts"]["fasta"],
        fai=config["ref_genome"] + config["file_exts"]["fai"],
        tumor_bam="{sampledir}/{tumor}" + config["file_exts"]["bam"],
        tumor_bai="{sampledir}/{tumor}" + config["file_exts"]["bai"],
        normal_bam="{sampledir}/{normal}" + config["file_exts"]["bam"],
        normal_bai="{sampledir}/{normal}" + config["file_exts"]["bai"]
    output:
        os.path.join("{sampledir}", get_outdir("lumpy"), \
            "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("lumpy")
    resources:
        mem_mb=get_maxmem("lumpy")
    shell:
        """
        echo {input} > {output}
        """
#    shell:
#        """
#        lumpyexpress -B {input.tumor_bam},{input.normal_bam} -o lumpy.vcf
#        """

rule lumpy:
    input:
        fasta=config["ref_genome"] + config["file_exts"]["fasta"],
        fai=config["ref_genome"] + config["file_exts"]["fai"],
        tumor_bam="{base_dir}/{tumor}" + config["file_exts"]["bam"],
        tumor_bai="{base_dir}/{tumor}" + config["file_exts"]["bai"],
        normal_bam="{base_dir}/{normal}" + config["file_exts"]["bam"],
        normal_bai="{base_dir}/{normal}" + config["file_exts"]["bai"]
    output:
        os.path.join("{base_dir}", get_outdir("lumpy"), \
            "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads: 1
    shell:
        """
        echo {input} > {output}
        """
#    shell:
#        """
#        lumpyexpress -B {input.tumor_bam},{input.normal_bam} -o lumpy.vcf
#        """


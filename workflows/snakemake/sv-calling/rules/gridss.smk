rule gridss:
    input:
        fasta=config["ref_genome"] + config["file_exts"]["fasta"],
        fai=config["ref_genome"] + config["file_exts"]["fai"],
        tumor_bam="{sampledir}/{tumor}" + config["file_exts"]["bam"],
        tumor_bai="{sampledir}/{tumor}" + config["file_exts"]["bai"],
        normal_bam="{sampledir}/{normal}" + config["file_exts"]["bam"],
        normal_bai="{sampledir}/{normal}" + config["file_exts"]["bai"]
        #TODO: Add BWA 0.6.x index files (.amb, .ann, .bwt, .pac, .sa)
    output:
        os.path.join("{sampledir}", get_outdir("gridss"), \
            "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("gridss")
    resources:
        mem_mb=get_maxmem("gridss")
    shell:
        """
        echo {input} > {output}
        """
#    shell:
#        """
#        gridss -Xmx31g gridss.CallVariants WORKER_THREADS={threads} \
#            REFERENCE_SEQUENCE={input.fasta} \
#            INPUT={input.normal_bam} \
#            INPUT={input.tumor_bam} \
#            OUTPUT=gridss.vcf \
#            ASSEMBLY=assembly.bam \
#            WORKING_DIR={output} \
#            TMP_DIR={outpu}
#        """

rule manta:
    input:
        fasta=config["ref_genome"] + config["file_exts"]["fasta"],
        fai=config["ref_genome"] + config["file_exts"]["fai"],
        tumor_bam="{sampledir}/{tumor}" + config["file_exts"]["bam"],
        tumor_bai="{sampledir}/{tumor}" + config["file_exts"]["bai"],
        normal_bam="{sampledir}/{normal}" + config["file_exts"]["bam"],
        normal_bai="{sampledir}/{normal}" + config["file_exts"]["bai"]
    output:
        os.path.join("{sampledir}", get_outdir("manta"), \
            "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("manta")
    resources:
        mem_mb=get_maxmem("manta")
    shell:
        """
        echo {input} > {output}
        """
#    shell:
#        """
#        configManta.py --runDir "{wildcards.sampledir}" \
#            --reference "{input.fasta}" \
#            --tumorBam "{input.tumor_bam}" \
#            --normalBam "{input.normal_bam}"
#        cd "{wildcards.sampledir}" && ./runWorkflow.py --quiet -m local -j {threads}
#        date "+%Y-%m-%d %H:%M:%S" > "{output}"
#        """

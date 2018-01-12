rule manta:
    input:
        fasta=get_fasta(),
        fai=get_faidx()[0],
        tumor_bam="{sampledir}/{tumor}" + get_filext("bam"),
        tumor_bai="{sampledir}/{tumor}" + get_filext("bam_idx"),
        normal_bam="{sampledir}/{normal}" + get_filext("bam"),
        normal_bai="{sampledir}/{normal}" + get_filext("bam_idx")
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

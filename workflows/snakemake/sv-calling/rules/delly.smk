rule delly:
    input:
        fasta=get_fasta(),
        fai=get_faidx()[0],
        tumor_bam="{sampledir}/{tumor}" + get_filext("bam"),
        tumor_bai="{sampledir}/{tumor}" + get_filext("bai"),
        normal_bam="{sampledir}/{normal}" + get_filext("bam"),
        normal_bai="{sampledir}/{normal}" + get_filext("bai")
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

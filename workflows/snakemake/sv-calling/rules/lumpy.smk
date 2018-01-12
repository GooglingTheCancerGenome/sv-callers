rule lumpy:
    input:
        fasta=get_fasta(),
        fai=get_faidx()[0],
        tumor_bam="{sampledir}/{tumor}" + get_filext("bam"),
        tumor_bai="{sampledir}/{tumor}" + get_filext("bam_idx"),
        normal_bam="{sampledir}/{normal}" + get_filext("bam"),
        normal_bai="{sampledir}/{normal}" + get_filext("bam_idx")
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

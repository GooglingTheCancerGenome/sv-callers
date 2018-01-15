rule delly_INV:
    input:
        fasta=get_fasta(),
        fai=get_faidx()[0],
        tumor_bam="{sampledir}/{tumor}" + get_filext("bam"),
        tumor_bai="{sampledir}/{tumor}" + get_filext("bam_idx"),
        normal_bam="{sampledir}/{normal}" + get_filext("bam"),
        normal_bai="{sampledir}/{normal}" + get_filext("bam_idx")
    params:
        outdir=os.path.join("{sampledir}", get_outdir("delly"))
    output:
        log=os.path.join("{sampledir}", get_outdir("delly"), \
            "{tumor}-{normal}_{sv_type, INV}.log")
    conda:
        "../environment.yaml"
    shell:
        delly_cmd()

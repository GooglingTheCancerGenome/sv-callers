rule delly:
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
            "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("delly")
    resources:
        mem_mb=get_maxmem("delly")
    shell:
        """
        if [ "{config[echo_run]}" = "1" ]; then
            echo "{input}" > "{output}"
        else
            # TODO: run all SV types in parallel
            delly call -t DUP \
                -g "{input.fasta}" \
                -o "{params}/delly-DUP.bcf" \
                "{input.tumor_bam}" "{input.normal_bam}" 2>&1
            # TODO:
            # delly filter?
            # bcf2vcf
            date "+%Y-%m-%d %H:%M:%S" > "{output}"
        fi
        """

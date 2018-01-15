rule manta:
    input:
        fasta=get_fasta(),
        fai=get_faidx()[0],
        tumor_bam="{sampledir}/{tumor}" + get_filext("bam"),
        tumor_bai="{sampledir}/{tumor}" + get_filext("bam_idx"),
        normal_bam="{sampledir}/{normal}" + get_filext("bam"),
        normal_bai="{sampledir}/{normal}" + get_filext("bam_idx")
    params:
        outdir=os.path.join("{sampledir}", get_outdir("manta"))
    output:
        log=os.path.join("{sampledir}", get_outdir("manta"), \
            "{tumor}-{normal}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("manta")
    resources:
        mem_mb=get_maxmem("manta")
    shell:
        """
        if [ "{config[echo_run]}" = "1" ]; then
            echo "{input}" > "{output}"
        else
            configManta.py \
                --runDir "{params}" \
                --reference "{input.fasta}" \
                --tumorBam "{input.tumor_bam}" \
                --normalBam "{input.normal_bam}" 2>&1
            cd "{params}" && \
            ./runWorkflow.py --quiet -m local -j {threads} 2>&1
            # TODO: rename/unpack outfile 'somaticSV.vcf.gz' to 'manta*.vcf'
            date "+%Y-%m-%d %H:%M:%S" > "{output}"
        fi
        """

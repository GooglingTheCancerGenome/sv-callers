rule manta:
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = "{path}/{tumor}" + get_filext("bam"),
        tumor_bai = "{path}/{tumor}" + get_filext("bam_idx"),
        normal_bam = "{path}/{normal}" + get_filext("bam"),
        normal_bai = "{path}/{normal}" + get_filext("bam_idx")
    output:
        log = os.path.join("{path}", "{tumor}--{normal}", get_outdir("manta"),
                           "{rule}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("manta")
    resources:
        mem_mb = get_memory("manta"),
        tmp_mb = get_tmpspace("manta")
    shell:
        """
        OUTDIR="$(dirname "{output}")"

        # run dummy or real job
        if [ "{config[echo_run]}" -eq "1" ]; then
            echo "{input}" > "{output}"
        else
            configManta.py \
                --runDir "${{OUTDIR}}" \
                --reference "{input.fasta}" \
                --tumorBam "{input.tumor_bam}" \
                --normalBam "{input.normal_bam}" &&
            cd "${{OUTDIR}}" &&
            ./runWorkflow.py \
                --quiet \
                -m local \
                -j {threads} &&
            zcat results/variants/somaticSV.vcf.gz > manta.vcf 2>&1
            date "+%Y-%m-%d %H:%M:%S" > "$(basename "{output}")"
        fi
        """

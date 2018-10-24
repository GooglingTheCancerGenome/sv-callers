rule manta_p:  # paired-samples analysis
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = get_bam("{path}/{tumor}"),
        tumor_bai = get_bai("{path}/{tumor}"),
        normal_bam = get_bam("{path}/{normal}"),
        normal_bai = get_bai("{path}/{normal}")
    # params:
    #     excl_opt = '-x "%s"' % get_bed("manta") if get_bed("manta") else ""
    output:
        os.path.join("{path}/{tumor}--{normal}", get_outdir("manta"),
                     "manta" + get_filext("vcf"))
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("manta")
    resources:
        mem_mb = get_memory("manta"),
        tmp_mb = get_tmpspace("manta")
    shell:
        """
        set -x

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
            # SV quality filtering
            bcftools filter \
                -O v `# uncompressed VCF format` \
                -o "$(basename "{output}")" \
                -i "FILTER == 'PASS'" \
                results/variants/somaticSV.vcf.gz
        fi
        """

rule manta_s:  # single-sample analysis: germline or tumor-only
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        bam = get_bam("{path}/{sample}"),
        bai = get_bai("{path}/{sample}")
    params:
        # excl_opt = '-x "%s"' % get_bed("manta") if get_bed("manta") else ""
        bam_opt = "--tumorBam" if is_tumor_only() else "--bam",
        outfile = "tumorSV.vcf.gz" if is_tumor_only() else "candidateSV.vcf.gz"
    output:
        os.path.join("{path}/{sample}", get_outdir("manta"), "manta" +
                     get_filext("vcf"))
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("manta")
    resources:
        mem_mb = get_memory("manta"),
        tmp_mb = get_tmpspace("manta")
    shell:
        """
        set -x

        OUTDIR="$(dirname "{output}")"

        # run dummy or real job
        if [ "{config[echo_run]}" -eq "1" ]; then
            echo "{input}" > "{output}"
        else
            configManta.py \
                --runDir "${{OUTDIR}}" \
                --reference "{input.fasta}" \
                {params.bam_opt} "{input.bam}" &&
            cd "${{OUTDIR}}" &&
            ./runWorkflow.py \
                --quiet \
                -m local \
                -j {threads} &&
            bcftools convert \
                -O v `# uncompressed VCF format` \
                -o "$(basename "{output}")" \
                results/variants/{params.outfile}
        fi
        """

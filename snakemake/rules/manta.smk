rule manta_s:  # somatic mode
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = "{path}/{tumor}" + get_filext("bam"),
        tumor_bai = "{path}/{tumor}" + get_filext("bam_idx"),
        normal_bam = "{path}/{normal}" + get_filext("bam"),
        normal_bai = "{path}/{normal}" + get_filext("bam_idx")
    #params:
    #    excl_opt = "-XXX " + get_bed("manta") if get_bed("manta") else ""
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
            bcftools convert \
                -O v `# uncompressed VCF format` \
                -o "$(basename "{output}")" \
                results/variants/somaticSV.vcf.gz
        fi
        """

rule manta_g:  # germline mode
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = "{path}/{tumor}" + get_filext("bam"),
        tumor_bai = "{path}/{tumor}" + get_filext("bam_idx")
    output:
        os.path.join("{path}/{tumor}", get_outdir("manta"), "manta" +
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
                --tumorBam "{input.tumor_bam}" &&
            cd "${{OUTDIR}}" &&
            ./runWorkflow.py \
                --quiet \
                -m local \
                -j {threads} &&
            bcftools convert \
                -O v `# uncompressed VCF format` \
                -o "$(basename "{output}")" \
                results/variants/tumorSV.vcf.gz
        fi
        """

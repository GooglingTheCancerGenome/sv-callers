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
            "{tumor}-{normal}_{sv_type}.log")
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
            PREFIX="{params}/delly-{wildcards.sv_type}"
            TSV="sample_pairs.tsv"
            # somatic SV calling
            delly call -t {wildcards.sv_type} \
                -g "{input.fasta}" \
                -o "${{PREFIX}}.bcf" \
                "{input.tumor_bam}" "{input.normal_bam}" && \
            # somatic pre-filtering
            printf "{wildcards.tumor}\ttumor\n{wildcards.normal}\tcontrol" \
                > "${{TSV}}" && \
            delly filter -f somatic -s "${{TSV}}" -o "${{PREFIX}}.pre.bcf" \
                "${{PREFIX}}.bcf" && \
            # BCF to VCF format conversion
            bcftools view "${{outfile_prefix}}.pre.bcf" -O v \
                -o "${{PREFIX}}.vcf" 2>&1
            date "+%Y-%m-%d %H:%M:%S" > "{output}"
        fi
        """

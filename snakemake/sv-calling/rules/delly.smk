rule delly:
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = "{sampledir}/{tumor}" + get_filext("bam"),
        tumor_bai = "{sampledir}/{tumor}" + get_filext("bam_idx"),
        normal_bam = "{sampledir}/{normal}" + get_filext("bam"),
        normal_bai = "{sampledir}/{normal}" + get_filext("bam_idx")
    output:
        log = os.path.join("{sampledir}", get_outdir("delly"),
                           "{tumor}-{normal}", "delly-{sv_type}.log")
    conda:
        "../environment.yaml"
    threads:
        get_nthreads("delly")
    resources:
        mem_mb = get_memory("delly"),
        tmp_mb = get_tmpspace("delly")
    shell:
        """
        OUTDIR="$(dirname "{output}")"

        function get_samp_id() {{
            echo "$(samtools view -H ${{1}} | \
                   perl -lne 'print ${{1}} if /\sLB:(\S+)/' | \
                   head -n 1)"
        }}

        if [ "{config[echo_run]}" = "1" ]; then
            echo "{input}" > "{output}"
        else
            # use OpenMP for threaded jobs
            export OMP_NUM_THREADS={threads}
            #export OMP_PROC_BIND=true
            #export OMP_PLACES=threads

            # somatic SV calling
            PREFIX="${{OUTDIR}}/delly-{wildcards.sv_type}"
            TSV="${{OUTDIR}}/sample_pairs.tsv"
            delly call \
                -t "{wildcards.sv_type}" \
                -g "{input.fasta}" \
                -o "${{PREFIX}}.bcf" \
                -q 1 `# min.paired-end mapping quality` \
                -s 9 `# insert size cutoff, DELs only` \
                "{input.tumor_bam}" "{input.normal_bam}" &&
            # somatic pre-filtering
            #   create sample list
            TID=$(get_samp_id "{input.tumor_bam}")
            CID=$(get_samp_id "{input.normal_bam}")
            printf "${{TID}}\ttumor\n${{CID}}\tcontrol\n" > ${{TSV}} &&
            delly filter \
                -f somatic \
                -s "${{TSV}}" \
                -o "${{PREFIX}}.pre.bcf" \
                "${{PREFIX}}.bcf" &&
            # BCF to VCF format conversion
            bcftools view \
                -O v `# VCF format` \
                -o "${{PREFIX}}.vcf" \
                "${{PREFIX}}.pre.bcf" 2>&1
            # TODO: merge SV VCF files
            date "+%Y-%m-%d %H:%M:%S" > "{output}"
        fi
        """

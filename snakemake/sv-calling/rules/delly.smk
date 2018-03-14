rule delly:
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = "{path}/{tumor}" + get_filext("bam"),
        tumor_bai = "{path}/{tumor}" + get_filext("bam_idx"),
        normal_bam = "{path}/{normal}" + get_filext("bam"),
        normal_bai = "{path}/{normal}" + get_filext("bam_idx")
    output:
        vcf = [os.path.join("{path}", "{tumor}--{normal}", get_outdir("delly"),
               "{rule}-" + sv + ".vcf")
               sv for sv in config["sv_callers"]["delly"]["sv_types"]]

    conda:
        "../environment.yaml"
    threads:
        get_nthreads("delly")
    resources:
        mem_mb = get_memory("delly"),
        tmp_mb = get_tmpspace("delly")
    shell:
        """
        set -x
        OUTDIR="$(dirname "{output}")"

        # fetch sample ID from a BAM file
        function get_samp_id() {{
            echo "$(samtools view -H ${{1}} | \
                   perl -lne 'print ${{1}} if /\sLB:(\S+)/' | \
                   head -n 1)"
        }}

        # run dummy or real job
        if [ "{config[echo_run]}" -eq "1" ]; then
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
                "${{PREFIX}}.pre.bcf"
            # TODO: merge SV VCF files
            #date "+%Y-%m-%d %H:%M:%S" > "{output}"
        fi
        """

rule delly_merge:
    input:
        [os.path.join("{path}", "{tumor}--{normal}", "{outdir}",
         "delly-" + sv + ".vcf")
         sv for sv in config["sv_callers"]["delly"]["sv_types"]]
    output:
        os.path.join("{path}", "{tumor}--{normal}", "{outdir}", "delly.vcf")
    shell:
        """
        set -x
        bcftools merge \
            -O v \
            -o "{output}" \
            "{input}"
        """

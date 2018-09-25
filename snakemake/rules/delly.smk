rule delly_p:  # paired-samples analysis
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        tumor_bam = get_bam("{path}/{tumor}"),
        tumor_bai = get_bai("{path}/{tumor}"),
        normal_bam = get_bam("{path}/{normal}"),
        normal_bai = get_bai("{path}/{normal}")
    params:
        excl_opt = '-x "%s"' % get_bed("delly") if get_bed("delly") else ""
    output:
        os.path.join("{path}/{tumor}--{normal}", get_outdir("delly"),
                     "delly-{sv_type}" + get_filext("bcf"))
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
        PREFIX="${{OUTDIR}}/$(basename "{output}" .bcf).unfiltered"
        TSV="${{OUTDIR}}/sample_pairs.tsv"

        # fetch sample ID from a BAM file
        function get_samp_id() {{
            echo "$(samtools view -H ${{1}} | \
                   perl -lne 'print ${{1}} if /\sSM:(\S+)/' | \
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

            # SV calling
            delly call \
                -t "{wildcards.sv_type}" \
                -g "{input.fasta}" \
                -o "${{PREFIX}}.bcf" \
                -q 1 `# min.paired-end mapping quality` \
                -s 9 `# insert size cutoff, DELs only` \
                {params.excl_opt} \
                "{input.tumor_bam}" "{input.normal_bam}" &&
            # somatic filtering
            #   create sample list
            TID=$(get_samp_id "{input.tumor_bam}")
            CID=$(get_samp_id "{input.normal_bam}")
            printf "${{TID}}\ttumor\n${{CID}}\tcontrol\n" > ${{TSV}} &&
            delly filter \
                -f somatic \
                -s "${{TSV}}" \
                -o "{output}" \
                "${{PREFIX}}.bcf"
        fi
        """

rule delly_s:  # single-sample analysis
    input:
        fasta = get_fasta(),
        fai = get_faidx()[0],
        bam = get_bam("{path}/{sample}"),
        bai = get_bai("{path}/{sample}")
    params:
        excl_opt = '-x "%s"' % get_bed("delly") if get_bed("delly") else ""
    output:
        os.path.join("{path}/{sample}", get_outdir("delly"), "delly-{sv_type}" +
                     get_filext("bcf"))
    conda:
        "../environment.yaml"
    threads: 1
    resources:
        mem_mb = get_memory("delly"),
        tmp_mb = get_tmpspace("delly")
    shell:
        """
        set -x
        OUTDIR="$(dirname "{output}")"

        # run dummy or real job
        if [ "{config[echo_run]}" -eq "1" ]; then
            echo "{input}" > "{output}"
        else
            # use OpenMP for threaded jobs
            export OMP_NUM_THREADS={threads}

            # SV calling
            delly call \
                -t "{wildcards.sv_type}" \
                -g "{input.fasta}" \
                -o "{output}" \
                -q 1 `# min.paired-end mapping quality` \
                -s 9 `# insert size cutoff, DELs only` \
                {params.excl_opt} \
                "{input.bam}"
        fi
        """

rule delly_merge:  # used by both modes
    input:
        [os.path.join("{path}/{tumor}--{normal}", get_outdir("delly"),
                      "delly-" + sv + get_filext("bcf"))
         for sv in config["callers"]["delly"]["sv_types"]]
        if config["mode"].startswith("p") is True else
        [os.path.join("{path}/{sample}", get_outdir("delly"), "delly-" + sv +
                      get_filext("bcf"))
         for sv in config["callers"]["delly"]["sv_types"]]
    output:
        os.path.join("{path}/{tumor}--{normal}", get_outdir("delly"), "delly" +
                     get_filext("vcf"))
        if config["mode"].startswith("p") is True else
        os.path.join("{path}/{sample}", get_outdir("delly"), "delly" +
                     get_filext("vcf"))
    conda:
        "../environment.yaml"
    threads: 1
    resources:
        mem_mb = 1024,
        tmp_mb = 0
    shell:
        """
        set -x

        # run dummy or real job
        if [ "{config[echo_run]}" -eq "1" ]; then
            cat {input} > "{output}"
        else
            # concatenate rather than merge BCF files
            bcftools concat \
               -a `# allow overlaps` \
               -O v `# uncompressed VCF format` \
               -o "{output}" \
               {input}
       fi
       """

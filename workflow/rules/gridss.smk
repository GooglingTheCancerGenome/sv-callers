rule gridss_p:  # paired-samples analysis
    input:
        fasta=get_fasta(),
        fai=get_faidx(),  # bwa index files also required
        tumor_bam=get_bam("{path}/{tumor}"),
        tumor_bai=get_bai("{path}/{tumor}"),
        normal_bam=get_bam("{path}/{normal}"),
        normal_bai=get_bai("{path}/{normal}"),
    params:
        excl_opt="BLACKLIST={}".format(get_bed()) if exclude_regions() else "",
    output:
        os.path.join(
            "{path}",
            "{tumor}--{normal}",
            get_outdir("gridss"),
            "gridss{}".format(config.file_exts.vcf),
        ),
    conda:
        "../envs/caller.yaml"
    threads: config.callers.gridss.threads
    resources:
        mem_mb=config.callers.gridss.memory,
        tmp_mb=config.callers.gridss.tmpspace,
    shell:
        """
        set -xe

        # if 'tmpspace' set to >0MB use TMPDIR otherwise use OUTDIR
        OUTDIR="$(dirname "{output}")"
        PREFIX="$(basename "{output}" .vcf)"
        OUTFILE="${{OUTDIR}}/${{PREFIX}}.unfiltered.vcf"
        TMP=$([ "{resources.tmp_mb}" -eq "0" ] &&
            echo "${{OUTDIR}}" || echo "${{TMPDIR}}")

        # set JVM max. heap size dynamically (in GB)
        # N.B. don't allocate >31G due to Compressed Oops and JDK-8029679
        MAX_HEAP=$(LC_ALL=C printf "%.f" $(bc <<< "scale=2; \
            {resources.mem_mb} / 1024 * .8")) # max. 80% of requested mem
        MAX_HEAP=$([ "${{MAX_HEAP}}" -gt "31" ] && echo "31g" ||
            echo "${{MAX_HEAP}}g")
        export _JAVA_OPTIONS="-Djava.io.tmpdir=${{TMP}} -Xmx${{MAX_HEAP}}"

        # run dummy or real job
        if [ "{config.echo_run}" -eq "1" ]; then
            echo "{input}" "${{TMP}}" > "{output}"
        else
            # clean-up outdir prior to SV calling
            rm -fr ${{OUTDIR}}/*gridss*
            gridss gridss.CallVariants \
                WORKER_THREADS={threads} \
                REFERENCE_SEQUENCE="{input.fasta}" \
                {params.excl_opt} \
                INPUT="{input.normal_bam}" \
                INPUT="{input.tumor_bam}" \
                OUTPUT="${{OUTFILE}}" \
                ASSEMBLY="${{OUTDIR}}/gridss_assembly.bam" \
                WORKING_DIR="${{TMP}}" \
                TMP_DIR="${{TMP}}/gridss.${{RANDOM}}"
            # somatic + SV quality filtering
            #   'normal' sample assumes index 0
            bcftools filter \
                -O v `# uncompressed VCF format` \
                -o "{output}" \
                -i "FORMAT/QUAL[0] == 0 && FILTER == '.'" \
                "${{OUTFILE}}"
        fi
        """


rule gridss_s:  # single-sample analysis
    input:
        fasta=get_fasta(),
        fai=get_faidx(),  # bwa index files also required
        bam=get_bam("{path}/{sample}"),
        bai=get_bai("{path}/{sample}"),
    params:
        excl_opt="BLACKLIST={}".format(get_bed()) if exclude_regions() else "",
    output:
        os.path.join(
            "{path}",
            "{sample}",
            get_outdir("gridss"),
            "gridss{}".format(config.file_exts.vcf),
        ),
    conda:
        "../envs/caller.yaml"
    threads: config.callers.gridss.threads
    resources:
        mem_mb=config.callers.gridss.memory,
        tmp_mb=config.callers.gridss.tmpspace,
    shell:
        """
        set -xe

        # if 'tmpspace' set to >0MB use TMPDIR otherwise use OUTDIR
        OUTDIR="$(dirname "{output}")"
        PREFIX="$(basename "{output}" .vcf)"
        OUTFILE="${{OUTDIR}}/${{PREFIX}}.unfiltered.vcf"
        TMP=$([ "{resources.tmp_mb}" -eq "0" ] &&
            echo "${{OUTDIR}}" || echo "${{TMPDIR}}")

        # set JVM max. heap size dynamically (in GB)
        # N.B. don't allocate >31G due to Compressed Oops and JDK-8029679
        MAX_HEAP=$(LC_ALL=C printf "%.f" $(bc <<< "scale=2; \
            {resources.mem_mb} / 1024 * .8")) # max. 80% of requested mem
        MAX_HEAP=$([ "${{MAX_HEAP}}" -gt "31" ] && echo "31g" ||
            echo "${{MAX_HEAP}}g")
        export _JAVA_OPTIONS="-Djava.io.tmpdir=${{TMP}} -Xmx${{MAX_HEAP}}"

        # run dummy or real job
        if [ "{config.echo_run}" -eq "1" ]; then
            echo "{input}" "${{TMP}}" > "{output}"
        else
            # clean-up outdir prior to SV calling
            rm -fr ${{OUTDIR}}/*gridss* &&
            gridss gridss.CallVariants \
                WORKER_THREADS={threads} \
                REFERENCE_SEQUENCE="{input.fasta}" \
                {params.excl_opt} \
                INPUT="{input.bam}" \
                OUTPUT="${{OUTFILE}}" \
                ASSEMBLY="${{OUTDIR}}/gridss_assembly.bam" \
                WORKING_DIR="${{TMP}}" \
                TMP_DIR="${{TMP}}/gridss.${{RANDOM}}"
            # SV quality filtering
            bcftools filter \
                -O v `# uncompressed VCF format` \
                -o "{output}" \
                -i "FILTER == '.'" \
                "${{OUTFILE}}"
        fi
        """

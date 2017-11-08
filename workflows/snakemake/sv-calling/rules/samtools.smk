rule samtools_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bai"
    threads: 8
    conda:
        "envs/sv_callers.yaml"
    shell:
        "samtools index -@ {threads} {input} {output}"

#
# This workflow executes samtools indexing of *.bam files
# provided in a config file.
#
# Example usage on SGE cluster:
#   snakemake --cluster "qsub -cwd -V -pe threaded {threads} -l h_rt=00:00:10 -l h_vmem=1G" \
#     --use-conda --jobs 2 --latency-wait 30 -s samtools.snakefile
#

configfile: "conf/samtools.yaml"

rule all:
    input:
        expand("{sample}.bai", sample=list(chain(*config["samples"])))

rule samtools_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bai"
    threads: 8
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} {output}"

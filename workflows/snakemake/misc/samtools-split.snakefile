#
# A simple workflow to split a whole-genome BAM file(s) into per-chromosome
# BAM files with proper read pairs only.
#
# Example usage on a SGE cluster:
#  snakemake --cluster "qsub -cwd -V -pe threaded {threads} -l h_rt=600 -l h_vmem=1G" \
#  --use-conda --jobs 10 -s samtools-split.snakefile {A,B}-{1,22}.bai
#
# This worklow will execute 10 jobs in parallel for two samples (A and B) and
# two chromosomes (1 and 2), and output four (indexed) BAM files.
#

rule samtools_view:
    input:
        bam="{sample}.bam",
        bai="{sample}.bai"
    output:
        "{sample}-{chr}.bam"
    threads: 8
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view -bh -f 2 -@ {threads} -o {output} {input.bam} {wildcards.chr}"

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

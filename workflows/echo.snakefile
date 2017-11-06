#
# A simple echo workflow using input parameters from a config file.
#
# Usage:
#   * on SGE cluster
#   snakemake --cluster "qsub -cwd -V -pe threaded {threads} -l h_rt=00:00:10 -l h_vmem=500M" \
#     --use-conda --jobs 2 --latency-wait 30 -s echo.snakefile
#
#   * on SLURM cluster
#   snakemake --cluster "sbatch --workdir=. --export=ALL -n {threads} -t 10 --mem=500" \
#     --use-conda --jobs 2 --latency-wait 30 -s echo.snakefile
#

configfile: "conf/echo.yaml"

rule all:
    input:
        expand("test/{sample}.bai", sample=config["samples"])

rule echo:
    input:
        bam="test/{sample}.bam"
    output:
        bai="test/{sample}.bai"
    threads: 1
    shell:
        "echo {input.bam} > {output.bai}"

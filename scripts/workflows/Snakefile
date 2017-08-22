# 
# A simple workflow to select proper reads mapped to a chromosome.
#
# To execute this workflow on a SGE cluster use:
#   snakemake --cluster "qsub -pe threaded {threads} -cwd -V" --jobs 1 samtools_index
#

CHR = 22
BAM_WG = "VCAP_dedup.realigned.bam"
BAI_WG = os.path.splitext(BAM_WG)[0] + ".bai"
BAM_CHR = "proper_pairs/chr%d.bam" % CHR
BAI_CHR = os.path.splitext(BAM_CHR)[0] + ".bai"

rule samtools_view:
    input:
        BAM_WG,
        BAI_WG
    output:
        BAM_CHR
    threads: 8
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view -bh -f 2 -@ {threads} -o {output} {input[0]} {CHR}"

rule samtools_index:
    input:
        BAM_CHR
    output:
        BAI_CHR
    threads: 8
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} {output}"

rule all:
    input:
        bai=expand("test/{sample}.bai", sample=config["samples"])

rule echo:
    input:
        bam="test/{sample}.bam"
    output:
        bai="test/{sample}.bai"
    shell:
        "echo {input.bam} {output.bai}"


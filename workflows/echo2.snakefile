#
# This workflow writes pairs of I/O files according samples in a config file.
#

configfile: "conf/echo2.yaml"

rule all:
    input:
        expand("{sample}.out", sample=config["samples"])

rule echo_pairs:
    input:
        "{id}_A.in",
        "{id}_B.in"
    output:
        "{id}_A.out",
        "{id}_B.out"
    shell:
        """
        echo A: {input[0]} > {output[0]}
        echo B: {input[1]} > {output[1]}
        """

rule write_files:
    output:
        "{id}_A.in",
        "{id}_B.in"
    shell:
        "touch {output}"

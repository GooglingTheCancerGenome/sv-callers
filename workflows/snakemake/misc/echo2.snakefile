#
# This workflow writes output file per sample pair given a config file.
#

import pandas as pd

shell.executable("/bin/bash")
configfile: "conf/echo2.yaml"

def get_samples(fname):
    return pd.read_csv(fname, sep='\t')

SAMPLES = get_samples(config['samples_file'])

rule all:
    input:
        [ "test/%s-%s.out" % (row.A, row.B) for row in SAMPLES.itertuples() ]

rule prep_input:
    output:
        expand("test/{sample}.in", sample=SAMPLES.values.flatten())
    shell:
        "touch {output}"

rule echo:
    input:
        A=expand("test/{A}.in", A=SAMPLES.A),
        B=expand("test/{B}.in", B=SAMPLES.B)
    output:
        "test/{A}-{B}.out"
    shell:
        "echo {input} > {output}"

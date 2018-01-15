def delly_cmd():
    """Return DELLY related shell commands.
    """
    return """
        if [ "{{config[echo_run]}}" = "1" ]; then
            echo "{{input}}" > "{{output}}"
        else
            delly call -t {{wildcards.sv_type}} \
                -g "{input.fasta}" \
                -o "{{params}}/delly-{{wildcards.sv_type}}.bcf" \
                "{{input.tumor_bam}}" "{{input.normal_bam}}" 2>&1
            # TODO:
            # 1. delly filter -f somatic -o *.pre.bcf -s samples.tsv *.bcf
            # 2. bcftools view *.pre.bcf -o *.pre.vcf -O v
            date "+%Y-%m-%d %H:%M:%S" > "{{output}}"
        fi
    """

rule delly:
    threads:
        get_nthreads("delly")
    resources:
        mem_mb=get_maxmem("delly")

include: "delly-BND.smk"
include: "delly-DEL.smk"
include: "delly-DUP.smk"
include: "delly-INS.smk"
include: "delly-INV.smk"

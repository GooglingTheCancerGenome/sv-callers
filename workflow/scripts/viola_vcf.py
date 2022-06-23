import os
import shutil
import viola

vcf_in = str(snakemake.input)
vcf_org = vcf_in + '.org'
vcf_out = str(snakemake.output)
caller = str(snakemake.wildcards.prefix)

if caller == 'gridss':
    os.rename(vcf_in, vcf_org)
    with open(vcf_in, 'w') as new:
        with open(vcf_org, 'r') as org:
            for line in org:
                # FIX INFO field: change PARID to MATEID
                new.write(line.replace('PARID', 'MATEID'))
try:
    sv = viola.read_vcf(vcf_in, variant_caller=caller).breakend2breakpoint()
    sv.to_vcf(vcf_out)
except Exception:
    shutil.copyfile(vcf_in, vcf_out)

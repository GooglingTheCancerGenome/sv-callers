# Filter GRIDSS output for somatic SVs

from pysam import VariantFile


vcf_infile = 'gridss.vcf'
vcf_outfile = 'gridss_somatic.vcf'

vcf_in = VariantFile(vcf_infile)
vcf_out = VariantFile(vcf_outfile, 'w', header=vcf_in.header)

for record in vcf_in:
    normal = record.samples.keys()[0]  # assume 1st sample is 'normal'
    if record.samples[normal]["QUAL"] == 0:
        vcf_out.write(record)

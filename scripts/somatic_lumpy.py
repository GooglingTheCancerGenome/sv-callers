# Outputs the somatic variants Lumpy found

from pysam import VariantFile

vcf_infile = 'lumpy.vcf'
vcf_outfile = 'lumpy.somatic.vcf'

vcf_in = VariantFile(vcf_infile)
vcf_out = VariantFile(vcf_outfile, 'w', header=vcf_in.header)

for record in vcf_in:
    if record.samples['COLO829T']['SU'] > 0 and record.samples['COLO829R']['SU'] == 0:
        vcf_out.write(record)

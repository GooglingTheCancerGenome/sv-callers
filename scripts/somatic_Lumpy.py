# Outputs the somatic variants Lumpy found

from pysam import VariantFile

vcf_input_file = wd + 'lumpy.vcf'
vcf_output_file = wd + 'lumpy.somatic.vcf'

vcf_in = VariantFile(vcf_input_file)
vcf_out = VariantFile(vcf_output_file, 'w', header=vcf_in.header)

for record in vcf_in:
    if record.samples['COLO829T']['SU']>0 and record.samples['COLO829R']['SU']==0:
        vcf_out.write(record)

# Outputs the somatic variants GRIDSS found

from pysam import VariantFile

vcf_input_file = 'somatic.sv.vcf'
vcf_output_file = 'somatic.vcf'

vcf_in = VariantFile(vcf_input_file)
vcf_out = VariantFile(vcf_output_file, 'w', header=vcf_in.header)

for record in vcf_in:
    #Normal is the first sample
    normal = record.samples.keys()[0]
    assert normal[0] == 'N'
    if record.samples[normal]['QUAL']==0:
        vcf_out.write(record)

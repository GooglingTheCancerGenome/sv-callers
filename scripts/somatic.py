from pysam import VariantFile

vcf_input_file = 'somatic.sv.vcf'
vcf_output_file = 'somatic.vcf'

vcf_in = VariantFile(vcf_input_file)
vcf_out = VariantFile(vcf_output_file, 'w', header=vcf_in.header)

for record in vcf_in:
    if record.samples['N1.bam']['QUAL']==0:
        vcf_out.write(record)

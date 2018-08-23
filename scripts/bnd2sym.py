import argparse
import re
from time import time
from pysam import VariantFile
import locations

"""
SURVIVOR merge does not consider the bracket notation (SVTYPE=BND) properly for SVTPEs.
This script converts the bracket notation into symbolic notation (<SVTYPE>) where possible.
"""


def bnd2sym(vcf_input_file, vcf_output_file, sv_caller):

    print('Considering SV caller: %s' % sv_caller)

    vcf_in = VariantFile(vcf_input_file)
    vcf_out = VariantFile(vcf_output_file, 'w', header=vcf_in.header)

    if sv_caller == 'Manta':
        for record in vcf_in:

            locations.setupREs()

            altstr = str(record.alts[0])
            resultBP_sym = re.match(locations.__symbolicRE__, altstr)
            resultBP_gen = re.match(locations.__genomicRE__, altstr)
            resultBP_bnd = re.match(locations.__bpRE__, altstr)
            #
            vcf_out.write(record)

    else:
        for record in vcf_in:

            locations.setupREs()
            # if resultBP_gen:
            #     ref_seq = record.ref
            #     alt_seq = resultBP_gen.group(1)
            #     if len(ref_seq) > len(alt_seq) and record.info['SVTYPE'] =='DEL':
            #         #print('DEL: Ref=%s | Alt=%s' % (ref_seq, alt_seq))
            #         #print('%s %s %d' % (record.info['SVTYPE'], record.chrom, record.start))
            #         record.alts = ('<DEL>',)
            #         #print(record.info['SVLEN'])
            #         record.stop = record.start + abs(record.info['SVLEN'][0])
            #         #vcf_out.write(record)
            #
            #     elif len(ref_seq) < len(alt_seq) and record.info['SVTYPE'] =='INS':
            #         #print('INS: Ref=%s | Alt=%s' % (ref_seq, alt_seq))
            #         #print('%s %s %d' % (record.info['SVTYPE'], record.chrom, record.start))
            #         record.alts = ('<INS>',)
            #         record.stop = record.start
            #         #vcf_out.write(record)
            # elif resultBP_bnd:
            #     ct, chr2, pos2, indellen = locations.get_bnd_info(record, resultBP_bnd)
            #     #print('CT:%s, chr2:%s, pos2:%s, indellen:%s' % (ct, chr2, pos2, indellen))
            #     if record.chrom == chr2 and record.start < pos2 and ct == '3to5':
            #         record.alts = ('<DEL>',)
            #         record.stop = record.start + abs(record.info['SVLEN'][0])
            #         #vcf_out.write(record)

            altstr = str(record.alts[0])
            resultBP_sym = re.match(locations.__symbolicRE__, altstr)
            resultBP_gen = re.match(locations.__genomicRE__, altstr)
            resultBP_bnd = re.match(locations.__bpRE__, altstr)

            if resultBP_bnd:
                ct, chr2, pos2, indellen = locations.get_bnd_info(record, resultBP_bnd)
                # print('CT:%s, chr2:%s, pos2:%s, indellen:%s' % (ct, chr2, pos2, indellen))
                if record.chrom == chr2:
                    if record.start < pos2 and ct == '3to5':
                        record.alts = ('<DEL>',)
                        record.stop = pos2
                        record.info['SVTYPE'] = 'DEL'
                        vcf_out.write(record)
                    elif ct == '5to5' or ct == '3to3':
                        record.alts = ('<INV>',)
                        record.stop = pos2
                        record.info['SVTYPE'] = 'INV'
                        vcf_out.write(record)
                    else:
                        record.alts = ('<DUP>',)
                        record.stop = pos2
                        record.info['SVTYPE'] = 'DUP'
                        vcf_out.write(record)
                else:
                    vcf_out.write(record)
            else:
                vcf_out.write(record)


def main():

    # Test input and output
    mills2011in = 'sv_benchmark/input.na12878/processed/lumpy-Mills2011_sorted.vcf'
    mills2011out = 'sv_benchmark/input.na12878/processed/lumpy-Mills2011_sorted.sym.vcf'

    parser = argparse.ArgumentParser(description='Filter VCF files from Delly, GRIDSS, Manta and Lumpy')
    parser.add_argument('-i', '--input', type=str, default=mills2011in,
                        help="Specify input file (VCF)")
    parser.add_argument('-o', '--output', type=str, default=mills2011out,
                        help="Specify output (VCF)")
    parser.add_argument('-c', '--caller', type=str, default='Delly',
                        help="Specify SV caller")

    args = parser.parse_args()

    t0 = time()

    bnd2sym(vcf_input_file=args.input, vcf_output_file=args.output, sv_caller=args.caller)

    print(time() - t0)


if __name__ == '__main__':
    main()

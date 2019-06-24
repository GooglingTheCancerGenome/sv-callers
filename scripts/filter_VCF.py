import argparse
import re
from time import time
from intervaltree import IntervalTree
from collections import defaultdict
from pysam import VariantFile

'''
This python script filters the VCF entries by quality and using the DAC Blacklist regions
'''

__bpRE__ = None
__symbolicRE__ = None
__genomicRE__ = None

# DAC Blacklisted regions
# https://www.encodeproject.org/annotations/ENCSR636HFF/
black_file = '/Users/lsantuari/Documents/Data/ENCODE/ENCFF001TDO.bed'

# The functions setupREs, stdchrom, locFromBkpt and get_bnd_info are derived/modified from
# mergevcf by Jonathan Dursi (Simpson Lab) GitHub: https://github.com/ljdursi/mergevcf
def setupREs():
    global __genomicRE__
    global __symbolicRE__
    global __bpRE__
    if __genomicRE__ is None or __symbolicRE__ is None or __bpRE__ is None:
        __genomicRE__ = re.compile(r'[ACGTNacgtn]+')
        __symbolicRE__ = re.compile(r'.*<([A-Z:]+)>.*')
        __bpRE__ = re.compile(r'([ACGTNactgn\.]*)([\[\]])([a-zA-Z0-9\.]+:\d+)([\[\]])([ACGTNacgtn\.]*)')


def stdchrom(chrom):
    if chrom[0] == 'c':
        return chrom[3:]
    else:
        return chrom


def locFromBkpt(ref, pre, delim1, pair, delim2, post):
    # extract strand/orientation, position, and (possibly) inserted string
    # length from record eg [9:1678896[N.
    # result is result from re.match with the explicit BP regexp

    chpos = pair.split(':')
    # print(chpos[0])
    chr2 = stdchrom(chpos[0])
    pos2 = int(chpos[1])
    assert delim1 == delim2  # '['..'[' or ']'...']'
    joinedAfter = True
    extendRight = True
    connectSeq = ""

    if len(pre) > 0:
        connectSeq = pre
        joinedAfter = True
        assert len(post) == 0
    elif len(post) > 0:
        connectSeq = post
        joinedAfter = False

    if delim1 == "]":
        extendRight = False
    else:
        extendRight = True

    indellen = len(connectSeq) - len(ref)

    if joinedAfter:
        if extendRight:
            ct = '3to5'
        else:
            ct = '3to3'
    else:
        if extendRight:
            ct = '5to5'
        else:
            ct = '5to3'

    return ct, chr2, pos2, indellen


def get_bnd_info(record, resultBP):
    assert resultBP.group(2) == resultBP.group(4)
    ct, chr2, pos2, indellen = locFromBkpt(str(record.ref), resultBP.group(1),
                                           resultBP.group(2), resultBP.group(3), resultBP.group(4),
                                           resultBP.group(5))
    return (ct, chr2, pos2, indellen)


def load_blacklist_as_interval_tree():

    trees = defaultdict(IntervalTree)

    for line in open(black_file):
        chrom, begin, end = line.rstrip().split("\t")[:3]
        chrom = chrom[3:]
        name = '_'.join((chrom, begin, end))
        trees[chrom][int(begin):int(end)] = name

    return trees


def is_blacklisted(treedict, record):

    # Initially written for Manta VCF output
    setupREs()
    altstr = str(record.alts[0])
    resultBP_sym = re.match(__symbolicRE__, altstr)
    resultBP_gen = re.match(__genomicRE__, altstr)
    resultBP_bnd = re.match(__bpRE__, altstr)

    # print(record)

    if resultBP_sym:
        # print('__symbolicRE__')
        if record.chrom not in treedict.keys():
            return False
        else:
            if 'CIPOS' in record.info and 'CIEND' in record.info:
                cipos1_search = treedict[record.chrom][record.pos + record.info['CIPOS'][0]]
                cipos2_search = treedict[record.chrom][record.pos + record.info['CIPOS'][1]]
                ciend1_search = treedict[record.chrom][record.stop + record.info['CIEND'][0]]
                ciend2_search = treedict[record.chrom][record.stop + record.info['CIEND'][1]]
                if len(cipos1_search) + len(cipos2_search) + len(ciend1_search) + len(ciend2_search) == 0:
                    return False
                else:
                    # print('Blacklisted:')
                    # print(record)
                    # print('cipos1_search = '+str(cipos1_search))
                    # print('cipos2_search = ' + str(cipos2_search))
                    # print('ciend1_search = ' + str(ciend1_search))
                    # print('ciend2_search = ' + str(ciend2_search))
                    return True
            else:
                tree_search_pos = treedict[record.chrom][record.pos]
                tree_search_end = treedict[record.chrom][record.stop]
                if len(tree_search_pos) == 0 and len(tree_search_end) == 0:
                    return False
                else:
                    # print('Blacklisted:')
                    # print(record)
                    # print('tree_search_pos = '+str(tree_search_pos))
                    # print('tree_search_end = ' + str(tree_search_end))
                    return True

    elif resultBP_gen:
        # print('__genomicRE__')
        if record.chrom not in treedict.keys():
            return False
        else:
            tree_search_pos = treedict[record.chrom][record.pos]
            tree_search_end = treedict[record.chrom][record.stop]
            if len(tree_search_pos) == 0 and len(tree_search_end) == 0:
                return False
            else:
                # print('Blacklisted:')
                # print(record)
                # print('tree_search_pos = ' + str(tree_search_pos))
                # print('tree_search_end = ' + str(tree_search_end))
                return True

    elif resultBP_bnd:
        # print('__bpRE__')
        ct, chr2, pos2, indellen = get_bnd_info(record, resultBP_bnd)
        # print('CT:%s, chr2:%s, pos2:%s, indellen:%s' % (ct, chr2, pos2, indellen))

        if record.chrom not in treedict.keys() and chr2 not in treedict.keys():
            return False

        if record.chrom in treedict.keys():
            if 'CIPOS' in record.info:
                cipos1_search = treedict[record.chrom][record.pos + record.info['CIPOS'][0]]
                cipos2_search = treedict[record.chrom][record.pos + record.info['CIPOS'][1]]
                if len(cipos1_search) + len(cipos2_search) != 0:
                    # print('Blacklisted:')
                    # print(record)
                    # print('cipos1_search = ' + str(cipos1_search))
                    # print('cipos2_search = ' + str(cipos2_search))
                    return True
            else:
                tree_search_pos = treedict[record.chrom][record.pos]
                if len(tree_search_pos) != 0:
                    # print('Blacklisted:')
                    # print(record)
                    # print('tree_search_pos = ' + str(tree_search_pos))
                    return True
        if chr2 in treedict.keys():
            tree_search_end = treedict[chr2][pos2]
            if len(tree_search_end) != 0:
                # print('Blacklisted:')
                # print(record)
                # print('tree_search_end = ' + str(tree_search_end))
                return True

        return False


def filter_vcf(vcf_input_file, vcf_output_file, sv_caller):

    print('SV caller: %s' % sv_caller)

    vcf_in = VariantFile(vcf_input_file)
    vcf_out = VariantFile(vcf_output_file, 'w', header=vcf_in.header)

    treedict = load_blacklist_as_interval_tree()

    n_pass = 0
    n_removed = 0
    n_total = 0
    n_written = 0

    if sv_caller == 'GRIDSS':

        min_qual = 500
        filter_out_set = {'LOW_QUAL', 'ASSEMBLY_TOO_FEW_READ', 'ASSEMBLY_TOO_SHORT', \
                          'LOW_BREAKPOINT_SUPPORT', 'NO_ASSEMBLY', 'REF', 'SMALL_EVENT'}

        for record in vcf_in:
            n_total += 1
            if record.qual >= min_qual:
                if len(record.filter) == 0 or len(set(record.filter) & filter_out_set) == 0:
                    n_pass += 1
                    if not is_blacklisted(treedict, record):
                        vcf_out.write(record)
                        n_written += 1
                    else:
                        n_removed += 1

        print('Total:%d records' % n_total)
        print('PASS:%d records' % n_pass)
        print('PASS and blacklisted:%d records' % n_removed)
        print('Written:%d records' % n_written)


    elif sv_caller == 'Delly':
        for record in vcf_in:
            n_total+=1
            if 'LowQual' not in list(record.filter):
                n_pass += 1
                if not is_blacklisted(treedict, record):
                    vcf_out.write(record)
                    n_written += 1
                else:
                    n_removed += 1

        print('Total:%d records' % n_total)
        print('PASS:%d records' % n_pass)
        print('PASS and blacklisted:%d records' % n_removed)
        print('Written:%d records' % n_written)


    elif sv_caller == 'Manta':
        n_pass = 0
        n_removed = 0
        n_total = 0
        for record in vcf_in:
            n_total+=1
            if 'PASS' in list(record.filter):
                n_pass += 1
                if not is_blacklisted(treedict, record):
                    vcf_out.write(record)
                    n_written += 1
                else:
                    n_removed+=1

    elif sv_caller == 'NanoSV':
        n_pass = 0
        n_removed = 0
        n_total = 0
        for record in vcf_in:
            n_total += 1
            if 'PASS' in list(record.filter):
                n_pass += 1
                if not is_blacklisted(treedict, record):
                    vcf_out.write(record)
                    n_written += 1
                else:
                    n_removed += 1

    elif sv_caller == 'DeepSV':
        n_pass = 0
        n_removed = 0
        n_total = 0
        for record in vcf_in:
            n_total += 1
            if 'PASS' in list(record.filter):
                n_pass += 1
                if not is_blacklisted(treedict, record):
                    vcf_out.write(record)
                    n_written += 1
                else:
                    n_removed += 1

        print('Total:%d records' % n_total)
        print('PASS:%d records' % n_pass)
        print('PASS and blacklisted:%d records' % n_removed)
        print('Written:%d records' % n_written)

    elif sv_caller == 'Lumpy':
        n_removed = 0
        n_total = 0
        for record in vcf_in:
            n_total += 1
            if not is_blacklisted(treedict, record):
                vcf_out.write(record)
                n_written += 1
            else:
                n_removed += 1

        print('Total:%d records' % n_total)
        print('Blacklisted:%d records' % n_removed)
        print('Written:%d records' % n_written)

        # if record.FILTER == None:
        #     filters.add(item)
        # elif isinstance(record.FILTER, list):
        #     for item in record.FILTER:
        #         filters.add(item)


def main():

    #context = 'trio/NA12878'
    context = '/Users/lsantuari/Documents/Data/germline/patients/Patient1/SV/Filtered/Filtering_stages/SURVIVOR_merge/with_oversampling/input/'

    # working directory
    #mantain = context+'/SV/Unfiltered/manta.vcf'
    #mantaout = context+'/SV/Filtered/manta.flt.vcf'

    deepsvin = context+'deepsv.DEL.vcf'
    deepsvout = context + 'deepsv.DEL.filtered.vcf'

    parser = argparse.ArgumentParser(description='Filter VCF files from Delly, GRIDSS, Manta and Lumpy')
    parser.add_argument('-i', '--input', type=str, default=deepsvin,
                        help="Specify input file (VCF)")
    parser.add_argument('-o', '--output', type=str, default=deepsvout,
                        help="Specify output (VCF)")
    parser.add_argument('-c', '--caller', type=str, default='DeepSV',
                        help="Specify SV caller")

    args = parser.parse_args()

    t0 = time()
    filter_vcf(vcf_input_file=args.input, vcf_output_file=args.output, sv_caller=args.caller)
    print(time() - t0)


if __name__ == '__main__':
    main()

import argparse
import re
from time import time
from intervaltree import IntervalTree
from collections import defaultdict
from pysam import VariantFile

"""
Auxilliary file adapted from mergevcf (https://github.com/ljdursi/mergevcf, Jonathan Dursi, Simpson Lab)
"""

__bpRE__ = None
__symbolicRE__ = None
__genomicRE__ = None


def setupREs():
    global __genomicRE__
    global __symbolicRE__
    global __bpRE__
    if __genomicRE__ is None or __symbolicRE__ is None or __bpRE__ is None:
        __genomicRE__ = re.compile(r'(^[ACGTNacgtn]+$)')
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
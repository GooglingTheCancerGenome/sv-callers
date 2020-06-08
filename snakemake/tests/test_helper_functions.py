"""Test suite for helper_functions.py."""
import pytest
import helper_functions as hf


test_filext = [
    ('fasta', '.fasta'),
    ('fasta_idx', ['.fasta.fai', '.fasta.bwt', '.fasta.amb', '.fasta.ann', '.fasta.pac', '.fasta.sa']),
    ('bam', '.bam'),
    ('bam_idx', '.bam.bai'),
    ('vcf', '.vcf'),
    ('bcf', '.bcf'),
    ('bed', '.bed')
]

test_modes = [hf.config.mode.SINGLE_SAMPLE, hf.config.mode.PAIRED_SAMPLE]

test_callers = [
    ('manta', 'manta_out', 24, 16384, 0),
    ('delly', 'delly_out', 2, 8192, 0),
    ('lumpy', 'lumpy_out', 1, 32768, 0),
    ('gridss', 'gridss_out', 24, 63488, 0),
]

test_inter_output = [
    ("s", [
        'data/bam/3/T3/manta_out/survivor/manta.vcf',
        'data/bam/3/T3/delly_out/survivor/delly.vcf',
        'data/bam/3/T3/lumpy_out/survivor/lumpy.vcf',
        'data/bam/3/T3/gridss_out/survivor/gridss.vcf']),
    ("p", [
        'data/bam/3/T3--N3/manta_out/survivor/manta.vcf',
        'data/bam/3/T3--N3/delly_out/survivor/delly.vcf',
        'data/bam/3/T3--N3/lumpy_out/survivor/lumpy.vcf',
        'data/bam/3/T3--N3/gridss_out/survivor/gridss.vcf'])
]

test_output = [
    ("s", ["data/bam/3/T3/all.vcf"]),
    ("p", ["data/bam/3/T3--N3/all.vcf"])
]

test_args = [
    ('filter', ['"data/ENCFF001TDO.bed"', -1, -1, 0, -1]),
    ('merge', ['all.txt', 100, 1, 0, 0, 0, 0, 'all.vcf'])
]


def test_get_fasta():
    result = hf.get_fasta()
    expected = 'data/fasta/chr22' + hf.config.file_exts.fasta
    assert result == expected

def test_get_faidx():
    result = hf.get_faidx()
    expected = ['data/fasta/chr22' + e for e in hf.config.file_exts.fasta_idx]
    assert result == expected

def test_get_bam():
    result = hf.get_bam('T3')
    expected = 'T3' + hf.config.file_exts.bam
    assert result == expected

def test_get_bai():
    result = hf.get_bai('T3')
    expected = 'T3' + hf.config.file_exts.bam_idx
    assert result == expected

def test_file_is_empty__emptyfile_exception():
    with pytest.raises(Exception):
        hf.file_is_empty('data/bam/1/T1.bam')

@pytest.mark.parametrize("caller,outdir,nthreads,memory,tmpspace", test_callers)
def test_get_outdir(caller, outdir, nthreads, memory, tmpspace):
    result = hf.get_outdir(caller)
    assert result == outdir

def test_exclude_regions():
    result = hf.exclude_regions()
    assert result in (0, 1)

def test_get_bed():
    result = hf.get_bed()
    expected = 'data/ENCFF001TDO.bed'
    assert result == expected

def test_is_tumor_only():
    result = hf.is_tumor_only()
    assert result in (0, 1)

@pytest.mark.parametrize("test_input,expected", test_args)
def test_survivor_args(test_input, expected):
    result = hf.survivor_args(test_input)
    assert result == expected

@pytest.fixture(params=test_modes)
def set_mode(request):
    return request.param

@pytest.mark.parametrize("test_input,expected", test_inter_output)
def test_make_output(set_mode, test_input, expected):
    hf.config.mode = set_mode
    result = hf.make_output()
    if hf.config.mode == test_input:
        assert set(result) == set(expected)

@pytest.mark.parametrize("test_input,expected", test_output)
def test_make_all(set_mode, test_input, expected):
    hf.config.mode = set_mode
    result = hf.make_all()
    if hf.config.mode == test_input:
        assert set(result) == set(expected)

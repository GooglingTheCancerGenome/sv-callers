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


test_callers = [
    ('manta', 'manta_out', 24, 16384, 0),
    ('delly', 'delly_out', 2, 8192, 0),
    ('lumpy', 'lumpy_out', 1, 32768, 0),
    ('gridss', 'gridss_out', 24, 63488, 0),
]


test_mode_inter_output = [
    ('p', [
        'data/bam/3/T3--N3/manta_out/survivor/manta.vcf', 'data/bam/3/T3--N3/delly_out/survivor/delly.vcf',
        'data/bam/3/T3--N3/lumpy_out/survivor/lumpy.vcf', 'data/bam/3/T3--N3/gridss_out/survivor/gridss.vcf']),
    ('s', [
        'data/bam/3/T3/manta_out/survivor/manta.vcf', 'data/bam/3/T3/delly_out/survivor/delly.vcf',
        'data/bam/3/T3/lumpy_out/survivor/lumpy.vcf', 'data/bam/3/T3/gridss_out/survivor/gridss.vcf'])
]


test_args = [
    ('filter', ['"data/ENCFF001TDO.bed"', -1, -1, 0, -1]),
    ('merge', ['all.txt', 100, 1, 0, 0, 0, 0, 'all.vcf'])
]


test_mode_output = [
    ('p', {'data/bam/3/T3--N3/all.vcf'}),
    ('s', {'data/bam/3/T3/all.vcf'})
]


def test_get_callers():
    result = hf.get_callers()
    expected = ['manta', 'delly', 'lumpy', 'gridss']
    assert result == expected


@pytest.mark.parametrize("test_input,expected", test_filext)
def test_get_filext(test_input, expected):
    result = hf.get_filext(test_input)
    assert result == expected


def test_get_filext__unknownextension_exception():
    with pytest.raises(Exception) as e_info:
        hf.get_filext('nobodyknowns')


def test_get_fasta():
    result = hf.get_fasta()
    expected = 'data/fasta/chr22' + hf.get_filext('fasta')
    assert result == expected


def test_get_faidx():
    result = hf.get_faidx()
    expected = ['data/fasta/chr22' + e for e in hf.get_filext('fasta_idx')]
    assert result == expected


def test_get_bam():
    result = hf.get_bam('T3')
    expected = 'T3' + hf.get_filext('bam')
    assert result == expected


def test_get_bai():
    result = hf.get_bai('T3')
    expected = 'T3' + hf.get_filext('bam_idx')
    assert result == expected


def test_file_is_empty__emptyfile_exception():
    with pytest.raises(Exception) as exc_info:
        hf.file_is_empty('data/bam/1/T1.bam')


@pytest.mark.parametrize("caller,outdir,nthreads,memory,tmpspace", test_callers)
def test_get_outdir(caller, outdir, nthreads, memory, tmpspace):
    result = hf.get_outdir(caller)
    assert result == outdir


@pytest.mark.parametrize("caller,outdir,nthreads,memory,tmpspace", test_callers)
def test_get_nthreads(caller, outdir, nthreads, memory, tmpspace):
    result = hf.get_nthreads(caller)
    assert result == nthreads


@pytest.mark.parametrize("caller,outdir,nthreads,memory,tmpspace", test_callers)
def test_get_memory(caller, outdir, nthreads, memory, tmpspace):
    result = hf.get_memory(caller)
    assert result == memory


@pytest.mark.parametrize("caller,outdir,nthreads,memory,tmpspace", test_callers)
def test_get_tmpspace(caller, outdir, nthreads, memory, tmpspace):
    result = hf.get_tmpspace(caller)
    assert result == tmpspace


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


@pytest.mark.parametrize("test_input,expected", test_mode_inter_output)
def test_make_output(test_input, expected):
    hf.config['mode'] = test_input
    result = hf.make_output()
    assert result == expected


@pytest.mark.parametrize("test_input,expected", test_mode_output)
def test_make_all(test_input, expected):
    hf.config['mode'] = test_input
    result = hf.make_all()
    assert result == expected

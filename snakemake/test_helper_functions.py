import pytest

import helper_functions as hf


def test_get_callers():
    result = hf.get_callers()

    expected = ['manta', 'delly', 'lumpy', 'gridss']
    assert result == expected


@pytest.mark.parametrize("test_input,expected", [
    ('fasta', '.fasta'),
    ('fasta_idx', ['.fasta.fai', '.fasta.bwt', '.fasta.amb', '.fasta.ann', '.fasta.pac', '.fasta.sa']),
    ('bam', '.bam'),
    ('bam_idx', '.bam.bai'),
    ('vcf', '.vcf'),
])
def test_get_filext(test_input, expected):
    result = hf.get_filext(test_input)

    assert result == expected


def test_get_filext__unknownextension_exception():
    with pytest.raises(Exception) as exc_info:
        hf.get_filext('nobodyknows')
        assert str(exc_info) == "Unknown input file format 'nobodyknows'."


def test_get_fasta():
    result = hf.get_fasta()
    expected = 'data/fasta/chr22.fasta'
    assert result == expected


def test_get_faidx():
    result = hf.get_faidx()
    expected = [
        'data/fasta/chr22.fasta.fai',
        'data/fasta/chr22.fasta.bwt',
        'data/fasta/chr22.fasta.amb',
        'data/fasta/chr22.fasta.ann',
        'data/fasta/chr22.fasta.pac',
        'data/fasta/chr22.fasta.sa',
    ]
    assert result == expected


test_callers = [
    ('manta', 'manta_out', 24, 16384, 0),
    ('delly', 'delly_out', 2, 8192, 0),
    ('lumpy', 'lumpy_out', 1, 32768, 0),
    ('gridss', 'gridss_out', 24, 63488, 0),
]


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


def test_make_output():
    result = hf.make_output()
    expected = [
        'data/bam/3/T3--N3/manta_out/survivor/manta.vcf', 'data/bam/3/T3--N3/delly_out/survivor/delly.vcf',
        'data/bam/3/T3--N3/lumpy_out/survivor/lumpy.vcf', 'data/bam/3/T3--N3/gridss_out/survivor/gridss.vcf'
    ]
    assert result == expected

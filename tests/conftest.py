"""Shared pytest fixtures for testing."""
import os
import tempfile
import shutil
from pathlib import Path

import pytest


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    temp_path = tempfile.mkdtemp()
    yield Path(temp_path)
    shutil.rmtree(temp_path)


@pytest.fixture
def sample_config():
    """Sample configuration for testing."""
    return {
        "genome": "/path/to/reference.fasta",
        "exclusion_list": "/path/to/excluded_regions.bed",
        "samples": "samples.csv",
        "exclude_regions": 1,
        "file_exts": {
            "fasta": ".fasta",
            "fasta_idx": [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"],
            "bed": ".bed",
            "bam": ".bam",
            "bam_idx": ".bai",
            "vcf": ".vcf"
        },
        "callers": {
            "manta": {
                "outdir": "manta_output",
                "tumor_only": 0
            },
            "delly": {
                "outdir": "delly_output"
            },
            "lumpy": {
                "outdir": "lumpy_output"
            },
            "gridss": {
                "outdir": "gridss_output"
            }
        },
        "postproc": {
            "survivor": {
                "outdir": "survivor_output",
                "filter": {
                    "min_size": 50,
                    "max_size": 10000,
                    "min_freq": 0.1,
                    "min_sup": 5
                },
                "merge": {
                    "infile": "input.vcf",
                    "max_dist": 1000,
                    "min_sup": 2,
                    "use_type": 1,
                    "use_strand": 1,
                    "use_size": 1,
                    "min_size": 30,
                    "outfile": "merged.vcf"
                }
            }
        },
        "enable_callers": ["manta", "delly", "lumpy"],
        "mode": "single_sample"
    }


@pytest.fixture
def sample_fasta_file(temp_dir):
    """Create a sample FASTA file for testing."""
    fasta_file = temp_dir / "reference.fasta"
    fasta_file.write_text(">chr1\nATCGATCGATCGATCG\n>chr2\nGCTAGCTAGCTAGCTA\n")
    return fasta_file


@pytest.fixture
def sample_bed_file(temp_dir):
    """Create a sample BED file for testing."""
    bed_file = temp_dir / "excluded_regions.bed"
    bed_file.write_text("chr1\t100\t200\nexcluded_region\nchr2\t500\t600\texcluded_region\n")
    return bed_file


@pytest.fixture
def sample_bam_file(temp_dir):
    """Create a sample BAM file for testing."""
    bam_file = temp_dir / "sample.bam"
    bam_file.write_text("fake bam content")  # Not real BAM format, just for testing file existence
    return bam_file


@pytest.fixture
def sample_bai_file(temp_dir):
    """Create a sample BAI file for testing."""
    bai_file = temp_dir / "sample.bai"
    bai_file.write_text("fake bai content")  # Not real BAI format, just for testing file existence
    return bai_file


@pytest.fixture
def sample_csv_file(temp_dir):
    """Create a sample CSV file for testing."""
    csv_file = temp_dir / "samples.csv"
    csv_file.write_text("PATH,SAMPLE1,SAMPLE2\n/data/bam,sample1.bam,sample2.bam\n")
    return csv_file


@pytest.fixture
def mock_config_object():
    """Mock config object with attributes for testing."""
    class MockConfig:
        def __init__(self):
            self.genome = "/path/to/reference.fasta"
            self.exclusion_list = "/path/to/excluded_regions.bed"
            self.samples = "samples.csv"
            self.exclude_regions = 1
            
            # File extensions
            self.file_exts = type('obj', (object,), {
                'fasta': '.fasta',
                'fasta_idx': ['.fai', '.amb', '.ann', '.bwt', '.pac', '.sa'],
                'bed': '.bed',
                'bam': '.bam',
                'bam_idx': '.bai',
                'vcf': '.vcf'
            })()
            
            # Callers configuration
            self.callers = type('obj', (object,), {
                'manta': type('obj', (object,), {
                    'outdir': 'manta_output',
                    'tumor_only': 0
                })(),
                'delly': type('obj', (object,), {'outdir': 'delly_output'})(),
                'lumpy': type('obj', (object,), {'outdir': 'lumpy_output'})(),
                'gridss': type('obj', (object,), {'outdir': 'gridss_output'})()
            })()
            
            # Post-processing configuration
            self.postproc = type('obj', (object,), {
                'survivor': type('obj', (object,), {
                    'outdir': 'survivor_output',
                    'filter': type('obj', (object,), {
                        'min_size': 50,
                        'max_size': 10000,
                        'min_freq': 0.1,
                        'min_sup': 5
                    })(),
                    'merge': type('obj', (object,), {
                        'infile': 'input.vcf',
                        'max_dist': 1000,
                        'min_sup': 2,
                        'use_type': 1,
                        'use_strand': 1,
                        'use_size': 1,
                        'min_size': 30,
                        'outfile': 'merged.vcf'
                    })()
                })()
            })()
            
            self.enable_callers = ["manta", "delly", "lumpy"]
            self.mode = type('obj', (object,), {
                'PAIRED_SAMPLE': 'paired_sample'
            })()
            self.mode.PAIRED_SAMPLE = 'paired_sample'
    
    return MockConfig()


@pytest.fixture
def empty_file(temp_dir):
    """Create an empty file for testing error conditions."""
    empty = temp_dir / "empty.txt"
    empty.touch()
    return empty


@pytest.fixture
def non_existent_file(temp_dir):
    """Return path to a non-existent file for testing error conditions."""
    return temp_dir / "does_not_exist.txt"
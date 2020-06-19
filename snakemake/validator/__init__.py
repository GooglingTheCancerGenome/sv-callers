"""YAML validator for analysis.yaml file."""
import enum
import logging

from typing import Dict, List
from ruamel import yaml

import yatiml


class Mode(enum.Enum):
    SINGLE_SAMPLE = 's'
    PAIRED_SAMPLE = 'p'

    @classmethod
    def _yatiml_recognize(cls, node: yatiml.UnknownNode) -> None:
        err_msg = ("Mode must be either 's' or 'p'")
        if node.yaml_node.value not in ('s', 'p'):
            raise yatiml.RecognitionError(err_msg)

    @classmethod
    def _yatiml_savorize(cls, node: yatiml.Node) -> None:
        yaml_to_py = {v._value_: v._name_ for v in cls.__members__.values()}
        if node.is_scalar(str):
            node.set_value(yaml_to_py.get(node.get_value()))


class FileExtension:
    """
    FileExtension class to hold the 'file_exts' node.
    """
    def __init__(self, fasta: str, fasta_idx: List[str], bam: str, bam_idx: str,
                 bed: str, bcf: str, vcf: str) -> None:
        self.fasta = fasta
        self.fasta_idx = fasta_idx
        self.bam = bam
        self.bam_idx = bam_idx
        self.bed = bed
        self.bcf = bcf
        self.vcf = vcf


class Resource:
    def __init__(self, threads: int, memory: int, tmpspace: int, outdir: str) -> None:
        self.threads = threads
        self.memory = memory
        self.tmpspace = tmpspace
        self.outdir = outdir


class Manta(Resource):
    def __init__(self, threads: int, memory: int, tmpspace: int, outdir: str,
                 tumor_only: int) -> None:
        super().__init__(threads, memory, tmpspace, outdir)
        self.tumor_only = tumor_only


class Delly(Resource):
    def __init__(self, threads: int, memory: int, tmpspace: int, outdir: str,
                 sv_types: List[str]) -> None:
        super().__init__(threads, memory, tmpspace, outdir)
        self.sv_types = sv_types


class Lumpy(Resource):
    def __init__(self, threads: int, memory: int, tmpspace: int,
                 outdir: str) -> None:
        super().__init__(threads, memory, tmpspace, outdir)


class Gridss(Resource):
    def __init__(self, threads: int, memory: int, tmpspace: int,
                 outdir: str) -> None:
        super().__init__(threads, memory, tmpspace, outdir)


class Caller:
    def __init__(self, manta: Manta, delly: Delly, lumpy: Lumpy,
                 gridss: Gridss) -> None:
        self.manta = manta
        self.delly = delly
        self.lumpy = lumpy
        self.gridss = gridss


class SurvivorFilter:
    def __init__(self, min_size: int, max_size: int, min_freq: int,
                 min_sup: int) -> None:
        self.min_size = min_size
        self.max_size = max_size
        self.min_freq = min_freq
        self.min_sup = min_sup


class SurvivorMerge:
    def __init__(self, infile: str, max_dist: int, min_sup: int, use_type: int,
                 use_strand: int, use_size: int, min_size: int, outfile: str) \
                 -> None:
        self.infile = infile
        self.max_dist = max_dist
        self.min_sup = min_sup
        self.use_type = use_type
        self.use_strand = use_strand
        self.use_size = use_size
        self.min_size = min_size
        self.outfile = outfile


class Survivor(Resource):
    def __init__(self, threads: int, memory: int, tmpspace: int, outdir: str,
                 filter: SurvivorFilter, merge: SurvivorMerge) -> None:
        super().__init__(threads, memory, tmpspace, outdir)
        self.filter = filter
        self.merge = merge


class PostProcess:
    def __init__(self, survivor: Survivor) -> None:
        self.survivor = survivor


class Analysis:
    """
    Analysis class to hold all nodes.
    """
    def __init__(self, echo_run: int, enable_callers: List[str], mode: Mode,
                 genome: str, exclusion_list: str, exclude_regions: int,
                 file_exts: FileExtension, samples: str, callers: Caller,
                 postproc: PostProcess) -> None:
        self.echo_run = echo_run
        self.enable_callers = enable_callers
        self.mode = mode
        self.genome = genome
        self.exclusion_list = exclusion_list
        self.exclude_regions = exclude_regions
        self.file_exts = file_exts
        self.samples = samples
        self.callers = callers
        self.postproc = postproc


# Create loader
class MyLoader(yatiml.Loader):
    """
    MyLoader class.
    """
    pass


def load_configfile(yaml_file: str) -> Dict:
    with open(yaml_file, 'r') as conf:
        return yaml.load(conf, MyLoader)


yatiml.logger.setLevel(logging.DEBUG)
yatiml.add_to_loader(MyLoader, [Mode, FileExtension, Resource, Manta, Delly,
                     Lumpy, Gridss, Caller, SurvivorFilter, SurvivorMerge,
                     Survivor, PostProcess, Analysis])
yatiml.set_document_type(MyLoader, Analysis)

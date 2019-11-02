## Test data for _sv-callers_ workflow

This distribution contains data analyzed by the _sv-callers_ workflow (v1.1.0).
It includes the human reference genomes in FASTA, excluded regions in BED(PE),
structural variants (SVs) in VCF, workflow config files and Jupyter Notebooks
(with R codes) to analyze SV callsets. Large BAM files were not included but
can be downloaded from the NCBI server [1] and the ENA database [2].

[1] [NA12878 sample](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam).

[2] [COLO829 tumor](https://identifiers.org/ena.embl:ERX2765496) with
matched [normal](https://identifiers.org/ena.embl:ERX2765495) sample.

### 1. Install dependencies.

```bash
sudo apt-get install libopenblas-dev  # on Debian-based Linux distros

conda env create -f example_1/environment.yaml  # or example_2
conda activate example_1                        #
R BATCH -e "IRkernel::installspec()"            # enable R kernel
```

### 2. Run analyses.

```bash
jupyter notebook
```

Note: Each example works only in its own conda env due to some dependency conflicts.

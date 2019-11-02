## Test data for _sv-callers_ workflow

This distribution includes data analyzed by the _sv-callers_ workflow (v1.1.0)
in the single-sample (germline) and paired-sample (somatic) modes:

- GRCh37 and b37 human reference genomes (`.fasta`)
- excluded genomic regions (`.bed(pe)`)
  - [ENCODE:ENCFF001TDO](http://identifiers.org/encode/ENCFF001TDO)
  - [Layer et al. (2014)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84#ref-CR28)
- structural variants (SVs) detected by the workflow (`.vcf`)
- SV truth sets (`.bed(pe)`)
  - [Personalis/1000 Genomes Project](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed)
  data by [Parikh et al. (2016)](https://doi.org/10.1186/s12864-016-2366-2)
  - [PacBio/Moleculo](https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2014-15-6-r84/MediaObjects/13059_2013_3363_MOESM4_ESM.zip) data by [Layer et al. (2014)](https://doi.org/10.1186/gb-2014-15-6-r84)
- workflow samples (`.csv`) and config files (`.yaml`)
- short-read alignments are not included due to large sizes but are
  freely available for download (`.bam`)
  - NA12878 [sample](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam)
  - COLO829 [tumor sample](https://identifiers.org/ena.embl:ERX2765496) with matched
    [normal sample](https://identifiers.org/ena.embl:ERX2765495)
- Jupyter Notebooks to analyze SV callsets (`.ipynb`)


### 1. Install Jupyter Notebook and dependencies.

```bash
sudo apt-get install libopenblas-dev  # on Debian-based Linux distros
conda env create -f example_1/environment.yaml  # or create &
conda activate example_1                        # activate example_2 env
R BATCH -e "IRkernel::installspec()"            # enable R kernel
```

### 2. Run analyses.

```bash
jupyter notebook
```

Note: Currently, each example notebook works only in its own environment due to
some dependency conflicts.

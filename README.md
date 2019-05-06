# sv-callers

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1217111.svg)](https://doi.org/10.5281/zenodo.1217111)
[![Build Status](https://travis-ci.org/GooglingTheCancerGenome/sv-callers.svg?branch=master)](https://travis-ci.org/GooglingTheCancerGenome/sv-callers)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/5377b44904f84f9380b9d778c49bdb9e)](https://app.codacy.com/app/arnikz/sv-callers?utm_source=github.com&utm_medium=referral&utm_content=GooglingTheCancerGenome/sv-callers&utm_campaign=Badge_Grade_Dashboard)

Structural variants (SVs) are an important class of genetic variation implicated in a wide array of genetic diseases. _sv-callers_ is a _Snakemake_-based workflow that combines several state-of-the-art tools for detecting SVs in whole genome sequencing (WGS) data. The workflow is easy to use and deploy on any Linux-based machine. In particular, the workflow supports automated software deployment, easy configuration and addition of new analysis tools as well as enables to scale from a single computer to different HPC clusters with minimal effort.

### Dependencies

- python  (>=3.6)
- [conda](https://conda.io/)  (>=4.5)
- [snakemake](https://snakemake.readthedocs.io/)  (>=4.8)
- [xenon-cli](https://github.com/NLeSC/xenon-cli)  (3.0.0)

These will be installed by the workflow itself:

- SV callers
  - [Manta](https://github.com/Illumina/manta)  (1.1.0)
  - [DELLY](https://github.com/dellytools/delly)  (0.7.7)
  - [LUMPY](https://github.com/arq5x/lumpy-sv)  (0.2.13)
  - [GRIDSS](https://github.com/PapenfussLab/gridss)  (1.3.4)

- Post-processing
  - [BCFtools](https://github.com/samtools/bcftools)  (1.9)
  - [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)  (1.0.6)

**1. Clone this repo.**

```bash
git clone https://github.com/GooglingTheCancerGenome/sv-callers.git
cd sv-callers/snakemake
```

**2. Install dependencies.**

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh  # download Miniconda installer (with Python 3)
bash miniconda.sh  # install Conda (accept defaults)
# add the following line to ~/.bashrc & source env
# export PATH="$HOME/miniconda3/bin:$PATH"
source ~/.bashrc
conda update -y conda  # update Conda
conda create -y -n wf && source activate wf  # create & activate new env
conda install -y -c bioconda snakemake
conda install -y -c nlesc xenon-cli  # optional but recommended;)
```

**3. Configure the workflow.**

- **config files**:
  - `analysis.yaml` - analysis-specific settings (e.g., workflow mode, I/O files, SV callers, post-processing or resources used etc.)
  - `environment.yaml` - software dependencies and versions

- **input files**:
  - example data in `sv-callers/snakemake/data` directory
  - reference genome in `.fasta` (incl. index files)
  - "exclusion" list of genomic regions in `.bed` (optional)
  - WGS samples in `.bam` (incl. index files)
  - list of (paired) samples in `samples.csv`

- **output files**:
  - (filtered) SVs per caller and merged calls in `.vcf` (incl. index files)

**4. Execute the workflow.**

```bash
# 'dry' run only checks I/O files
snakemake -np

# 'vanilla' run (default) mimics the execution of SV callers by writing (dummy) VCF files
snakemake -C echo_run=1

```

Note: One sample or a tumor/normal pair generates eight SV calling jobs (i.e., 1 x Manta, 1 x LUMPY, 1 x GRIDSS and 5 x DELLY) and one post-processing job to merge DELLY call sets (per SV type) into one output file in [VCF](https://samtools.github.io/hts-specs/). See the workflow instance of [single-sample](https://github.com/GooglingTheCancerGenome/sv-callers/blob/dev/doc/sv-callers_single.svg) (germline) or [paired-sample](https://github.com/GooglingTheCancerGenome/sv-callers/blob/dev/doc/sv-callers_paired.svg) (somatic) analysis.

_Submit jobs to Grid Engine-based cluster_

```bash
snakemake -C echo_run=1 mode=p enable_callers="['manta','delly','lumpy','gridss']" --use-conda --latency-wait 30 --jobs 9 \
--cluster 'xenon scheduler gridengine --location local:// submit --name smk.{rule} --inherit-env --procs-per-node {threads} --start-single-process --max-run-time 1 --max-memory {resources.mem_mb} --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log' &>smk.log&
```

_Submit jobs to Slurm-based cluster_

```bash
snakemake -C echo_run=1 mode=p enable_callers="['manta','delly','lumpy','gridss']" --use-conda --latency-wait 30 --jobs 9 \
--cluster 'xenon scheduler slurm --location local:// submit --name smk.{rule} --inherit-env --procs-per-node {threads} --start-single-process --max-run-time 1 --max-memory {resources.mem_mb} --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log' &>smk.log&
```

To perform SV calling:
- overwrite (default) parameters directly in `analysis.yaml` or via the _snakemake_ CLI (use the `-C` argument)
  - set `echo_run=0`
  - choose between two workflow `mode`s: single- (`s`) or paired-sample (`p` default)
  - select one or more callers using `enable_callers` (default all: `"['manta','delly,'lumpy','gridss']"`)

- use `xenon` CLI to set:
  - `--max-run-time` of workflow jobs (in minutes)
  - `--temp-space` (optional, in MB)

- adjust compute requirements per SV caller according to the system used:
  - the number of `threads`, 
  - the amount of `memory`(in MB) or
  - the amount of temporary disk space or `tmpspace` (path in `TMPDIR` env variable) can be used for intermediate files by LUMPY and GRIDSS only.

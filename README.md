# sv-callers

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1217111.svg)](https://doi.org/10.5281/zenodo.1217111)
[![Published in PeerJ](https://img.shields.io/badge/published%20in-PeerJ-blue.svg)](https://doi.org/10.7717/peerj.8214)
[![Build Status](https://travis-ci.org/GooglingTheCancerGenome/sv-callers.svg?branch=dev)](https://travis-ci.org/GooglingTheCancerGenome/sv-callers)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/5377b44904f84f9380b9d778c49bdb9e)](https://app.codacy.com/app/arnikz/sv-callers?utm_source=github.com&utm_medium=referral&utm_content=GooglingTheCancerGenome/sv-callers&utm_campaign=Badge_Grade_Dashboard)

Structural variants (SVs) are an important class of genetic variation implicated in a wide array of genetic diseases. _sv-callers_ is a _Snakemake_-based workflow that combines several state-of-the-art tools for detecting SVs in whole genome sequencing (WGS) data. The workflow is easy to use and deploy on any Linux-based machine. In particular, the workflow supports automated software deployment, easy configuration and addition of new analysis tools as well as enables to scale from a single computer to different HPC clusters with minimal effort.

## Dependencies

-   [Python 3](https://www.python.org/)
-   [Conda](https://conda.io/) - package/environment management system
-   [Snakemake](https://snakemake.readthedocs.io/) - workflow management system
-   [Xenon CLI](https://github.com/NLeSC/xenon-cli) - command-line interface to compute and storage resources
-   [jq](https://stedolan.github.io/jq/) - command-line JSON processor (optional)
-   [YAtiML](https://github.com/yatiml/yatiml) - library for YAML type inference and schema validation

The workflow includes the following bioinformatics tools:

-   SV callers
    -   [Manta](https://github.com/Illumina/manta) (1.1.0)
    -   [DELLY](https://github.com/dellytools/delly) (0.7.7)
    -   [LUMPY](https://github.com/arq5x/lumpy-sv) (0.2.13)
    -   [GRIDSS](https://github.com/PapenfussLab/gridss) (1.3.4)

-   Post-processing
    -   [BCFtools](https://github.com/samtools/bcftools) (1.9)
    -   [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) (1.0.6)

The software dependencies and versions can be found in the conda `environment.yaml` files ([1](/environment.yaml), [2](/snakemake/environment.yaml)).

**1. Clone this repo.**

```bash
git clone https://github.com/GooglingTheCancerGenome/sv-callers.git
cd sv-callers
```

**2. Install dependencies.**

```bash
# download Miniconda3 installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# install Conda (respond by 'yes')
bash miniconda.sh
# update Conda
conda update -y conda
# create & activate new env with installed deps
conda env create -n wf -f environment.yaml
conda activate wf
cd snakemake
```

**3. Configure the workflow.**

-   **config files**:
    -   [`analysis.yaml`](/snakemake/analysis.yaml) - analysis-specific settings (e.g., workflow mode, I/O files, SV callers, post-processing or resources used etc.)
    -   [`environment.yaml`](/snakemake/environment.yaml) - software dependencies and versions

-   **input files**:
    -   example data in `sv-callers/snakemake/data` directory
    -   reference genome in `.fasta` (incl. index files)
    -   excluded regions in `.bed` (optional)
    -   WGS samples in `.bam` (incl. index files)
    -   list of (paired) samples in `samples.csv`

-   **output files**:
    -   (filtered) SVs per caller and merged calls in `.vcf` (incl. index files)

**4. Execute the workflow.**

_Locally_

```bash
# 'dry' run only checks I/O files
snakemake -np

# 'vanilla' run mimics the execution of SV callers by writing (dummy) VCF files
# if echo_run set to 1 (default) in analysis.yaml 
snakemake --cores
```

_Submit jobs to Slurm or GridEngine cluster_

```bash
SCH=slurm   # or gridengine
snakemake  --use-conda --latency-wait 30 --jobs 14 \
--cluster "xenon scheduler $SCH --location local:// submit --name smk.{rule} --inherit-env --cores-per-task {threads} --max-run-time 1 --max-memory {resources.mem_mb} --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log" &>smk.log&
```

Note: One sample or a tumor/normal pair generates eight SV calling jobs (i.e., 1 x Manta, 1 x LUMPY, 1 x GRIDSS and 5 x DELLY) and six post-processing jobs. See the workflow instance of [single-sample](doc/sv-callers_single.svg) (germline) or [paired-sample](doc/sv-callers_paired.svg) (somatic) analysis.

To perform SV calling:
-   edit (default) parameters in `analysis.yaml`
    -   set `echo_run` to `0`
    -   choose between two workflow `mode`s: single- (`s`) or paired-sample (`p` - default)
    -   select one or more callers using `enable_callers` (default all)

-   use `xenon` CLI to set:
    -   `--max-run-time` of workflow jobs (in minutes)
    -   `--temp-space` (optional, in MB)

-   adjust compute requirements per SV caller according to the system used:
    -   the number of `threads`, 
    -   the amount of `memory`(in MB),
    -   the amount of temporary disk space or `tmpspace` (path in `TMPDIR` env variable) can be used for intermediate files by LUMPY and GRIDSS only.

_Query job accounting information_

```bash
SCH=slurm   # or gridengine
xenon --json scheduler $SCH --location local:// list --identifier [jobID] | jq ...
```

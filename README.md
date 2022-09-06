# sv-callers

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1217111.svg)](https://doi.org/10.5281/zenodo.1217111)
[![Published in PeerJ](https://img.shields.io/badge/published%20in-PeerJ-blue.svg)](https://doi.org/10.7717/peerj.8214)
[![CI](https://github.com/GooglingTheCancerGenome/sv-callers/actions/workflows/ci.yaml/badge.svg?branch=master)](https://github.com/GooglingTheCancerGenome/sv-callers/actions/workflows/ci.yaml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/eaa33d7d090048898c112a4a87815479)](https://www.codacy.com/gh/GooglingTheCancerGenome/sv-callers/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=GooglingTheCancerGenome/sv-callers&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/eaa33d7d090048898c112a4a87815479)](https://www.codacy.com/gh/GooglingTheCancerGenome/sv-callers/dashboard?utm_source=github.com&utm_medium=referral&utm_content=GooglingTheCancerGenome/sv-callers&utm_campaign=Badge_Coverage)

Structural variants (SVs) are an important class of genetic variation implicated in a wide array of genetic diseases. _sv-callers_ is a _Snakemake_-based workflow that combines several state-of-the-art tools for detecting SVs in whole genome sequencing (WGS) data. The workflow is easy to use and deploy on any Linux-based machine. In particular, the workflow supports automated software deployment, easy configuration and addition of new analysis tools as well as enables to scale from a single computer to different HPC clusters with minimal effort.

## Dependencies

-   [Python](https://www.python.org/)
-   [Conda](https://conda.io/) - package/environment management system
-   [Snakemake](https://snakemake.readthedocs.io/) - workflow management system
-   [Xenon CLI](https://github.com/NLeSC/xenon-cli) - command-line interface to compute and storage resources
-   [jq](https://stedolan.github.io/jq/) - command-line JSON processor (optional)
-   [YAtiML](https://github.com/yatiml/yatiml) - library for YAML type inference and schema validation

The workflow includes the following bioinformatics tools:

-   SV callers
    -   [Manta](https://github.com/Illumina/manta)
    -   [DELLY](https://github.com/dellytools/delly)
    -   [LUMPY](https://github.com/arq5x/lumpy-sv)
    -   [GRIDSS](https://github.com/PapenfussLab/gridss)

-   Post-processing
    -   [BCFtools](https://github.com/samtools/bcftools)
    -   [Viola-SV](https://github.com/dermasugita/Viola-SV)
    -   [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)

The software dependencies can be found in the conda environment files: [[1]](/environment.yaml),[[2]](/workflow/envs/caller.yaml),[[3]](/workflow/envs/postproc.yaml).

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
# install Mamba
conda install -n base -c conda-forge -y mamba
# create a new environment with dependencies & activate it
mamba env create -n wf -f environment.yaml
conda activate wf
```

**3. Configure the workflow.**

-   **config files**:
    -   [`analysis.yaml`](/config/analysis.yaml) - analysis-specific settings (e.g., workflow mode, I/O files, SV callers, post-processing or resources used etc.)
    -   [`samples.csv`](/config/samples.csv) - list of (paired) samples

-   **input files**:
    -   example data in `workflow/data` directory
    -   reference genome in `.fasta` (incl. index files)
    -   excluded regions in `.bed` (optional)
    -   WGS samples in `.bam` (incl. index files)

-   **output files**:
    -   (filtered) SVs per caller and merged calls in `.vcf` (incl. index files)

**4. Execute the workflow.**

```bash
cd workflow
```

_Locally_

```bash
# 'dry' run only checks I/O files
snakemake -np

# 'vanilla' run if echo_run set to 1 (default) in analysis.yaml,
# it merely mimics the execution of SV callers by writing (dummy) VCF files;
# SV calling if echo_run set to 0
snakemake --use-conda --jobs

```

_Submit jobs to Slurm or GridEngine cluster_

```bash
SCH=slurm   # or gridengine
snakemake  --use-conda --latency-wait 30 --jobs \
--cluster "xenon scheduler $SCH --location local:// submit --name smk.{rule} --inherit-env --cores-per-task {threads} --max-run-time 1 --max-memory {resources.mem_mb} --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log" &>smk.log&
```

Note: One sample or a tumor/normal pair generates in total 18 SV calling and post-processing jobs. See the workflow instance of [single-sample](doc/sv-callers_single.svg) (germline) or [paired-sample](doc/sv-callers_paired.svg) (somatic) analysis.

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

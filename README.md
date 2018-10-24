# sv-callers

[![Build Status](https://travis-ci.org/GooglingTheCancerGenome/sv-callers.svg?branch=dev)](https://travis-ci.org/GooglingTheCancerGenome/sv-callers)

Structural variants (SVs) are an important class of genetic variation implicated in a wide array of genetic diseases. _sv-callers_ is a _Snakemake_-based workflow that combines several state-of-the-art tools for detecting SVs in whole genome sequencing (WGS) data. The workflow is easy to use and deploy on any Linux-based machine. In particular, the workflow supports automated software deployment, easy configuration and addition of new analysis tools as well as enables to scale from a single computer to different HPC clusters with minimal effort.

### Dependencies

- python (>=3.6)
- [conda](https://conda.io/) (>=4.5)
- [snakemake](https://snakemake.readthedocs.io/) (>=4.7)
- [xenon-cli](https://github.com/NLeSC/xenon-cli) (2.4)
- SV callers (installed via the [bioconda](https://bioconda.github.io/) channel):
  - [Manta](https://github.com/Illumina/manta) (1.1.0)
  - [DELLY](https://github.com/dellytools/delly) (0.7.7)
  - [LUMPY](https://github.com/arq5x/lumpy-sv) (0.2.13)
  - [GRIDSS](https://github.com/PapenfussLab/gridss) (1.3.4)

**1. Clone this repo.**

```bash
git clone https://github.com/GooglingTheCancerGenome/sv-callers.git
cd sv-callers/snakemake
```

**2. Install dependencies.**

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh  # python3
bash miniconda.sh  # install & add conda to your PATH
source ~/.bashrc
conda update -y conda  # update conda
conda create -y -n wf && source activate wf  # create & activate the environment
conda install -y -c bioconda snakemake
conda install -y -c nlesc xenon-cli  # optional but recommended;)
```

**3. Configure and execute the workflow.**

- **config files**:
  - `analysis.yaml` - analysis-specific settings (e.g., workflow mode, SV callers used etc.)
  - `environment.yaml` - software dependencies and versions
- **input files**:
  - example data in the `sv-callers/data` directory
  - reference genome in `.fasta` (incl. index files)
  - (paired) samples in `*.bam` (incl. index files)
  - list of (paired) samples for analysis in `samples.csv`
- **output files**: SVs in `.vcf` (incl. index files)

Note: One sample or tumor/normal sample pair generates eight SV calling jobs (i.e. 1 x Manta, 1 x LUMPY, 1 x GRIDSS and 5 x DELLY) and one post-processing job to merge DELLY call sets (per SV type) into one VCF file. A workflow instance can be found [here](https://github.com/GooglingTheCancerGenome/sv-callers/blob/master/doc/sv_calling_workflow.png).

```bash
# the 'dry' run only checks I/O files
snakemake -np
# the 'vanilla' run (default) doesn't execute the SV callers per se but rather writes the I/O file names into (dummy) VCF files
snakemake -C echo_run=1

```

_Submit jobs to Grid Engine-based cluster_

```bash
snakemake -C echo_run=1 mode=p enable_callers="['manta','delly','lumpy','gridss']" --latency-wait 30 --jobs 9 \
--cluster 'xenon scheduler gridengine --location local:// submit --name smk.{rule} --inherit-env --option parallel.environment=threaded --option parallel.slots={threads} --max-run-time 1 --max-memory {resources.mem_mb} --working-directory . --stderr stderr-\\\$JOB_ID.log --stdout stdout-\\\$JOB_ID.log' &>smk.log&
```

_Submit jobs to Slurm-based cluster_

```bash
snakemake -C echo_run=1 mode=p enable_callers="['manta','delly','lumpy','gridss']" --latency-wait 30 --jobs 9 \
--cluster 'xenon scheduler slurm --location local:// submit --name smk.{rule} --inherit-env --procs-per-node {threads} --start-single-process --max-run-time 1 --max-memory {resources.mem_mb} --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log' &>smk.log&
```

Settings for SV calling:
- `echo_run=0`
- set `--max-run-time` to a reasonable value (in minutes)
- add `--use-conda` to run jobs in a conda environment with pre-install callers
- select one or more callers using `enable_callers` (default all: `"['manta','delly,'lumpy','gridss']"`)
- choose between two workflow modes: `s` - single-sample or `p` - paired-samples analysis (default: `p`)

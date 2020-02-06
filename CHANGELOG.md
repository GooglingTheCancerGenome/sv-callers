# Changelog

## [1.1.1] - 2020-02-06
- added `environment.yaml` to install initial dependencies via Conda
  - updated Xenon CLI (v3.0.4): `--cores-per-task {threads}` replaces `--procs-per-node {threads}` and `--start-single-process` args
  - included jq JSON processor to query job accounting output by Xenon CLI (optional)
- updated Docker images with batch schedulers for CI testing
- fixed file error handling and unit tests

## [1.1.0] - 2019-05-06
- enabled both single-sample (germline) and paired-sample (somatic) SV analyses
- added post-processing steps (i.e., SV quality filtering and/or excluding regions with problematic mappability)
- updated Xenon CLI (v3.0.0) and included optional use of temporary disk space
- added Docker images with batch schedulers for CI testing
- bug fixes

## [1.0.0] - 2018-04-11
- initial release supporting SV detection in paired tumor/normal samples

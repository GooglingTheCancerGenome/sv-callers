# Changelog

## [1.2.2] - 2022-11-23
- hotfix (#78)

## [1.2.1] - 2022-09-06
- fixed named wildcards in s/p modes (#76)
- updated to Viola-SV v1.0.2

## [1.2.0] - 2022-06-23
- integrate Viola-SV into the post-processing step to harmonize SV types (#75)

## [1.1.3] - 2022-02-08
- fixed filtering of somatic SV types in DELLY (#68)
- updated testing infra
  - Travis CI -> GitHub Actions (#55)
  - Docker Hub -> GitHub Container Registry (#71)

## [1.1.2] - 2020-10-26
- update Xenon CLI (v3.0.5): fixed Slurm 19 issue (#40)
- speed up conda install (#49)

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

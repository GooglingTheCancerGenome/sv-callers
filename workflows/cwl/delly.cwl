#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "DELLY structural variant caller"
baseCommand: delly
arguments: ["call"]

hints:
  SoftwareRequirement:
    packages:
    - package: delly
      version:
      - "0.7.7"
      specs:
      - https://anaconda.org/bioconda/delly

inputs:
  - id: sv_type
    type: string
    inputBinding:
      position: 1
      prefix: "-t"
      separate: true
  - id: excluded_regions
    type: File
    inputBinding:
      position: 2
      prefix: "-x"
      separate: true
  - id: genome
    type: File
    inputBinding:
      position: 3
      prefix: "-g"
      separate: true
  - id: variants
    type: string
    inputBinding:
      position: 4
      prefix: "-o"
      separate: true
  - id: tumor_sample
    type: File
    inputBinding:
      position: 5
  - id: normal_sample
    type: File
    inputBinding:
      position: 6

outputs:
  - id: bcf_output
    type: File
    outputBinding:
      glob: "*.bcf"

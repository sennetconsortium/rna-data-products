cwlVersion: v1.2
class: CommandLineTool
label: Identify and score scenescent cells
requirements:
  DockerRequirement:
    dockerPull: sennet/rna-data-products-python:latest
baseCommand: /opt/deepscence.py

inputs:
  raw_h5ad_file:
    type: File
    inputBinding:
      position: 0
outputs:
  raw_h5ad_with_ds:
    type: File
    outputBinding:
      glob: expr.h5ad
  deepscence_plot:
    type: File
    outputBinding:
      glob: deepscence.pdf
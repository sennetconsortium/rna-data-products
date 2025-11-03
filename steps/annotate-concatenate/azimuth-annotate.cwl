cwlVersion: v1.0
class: CommandLineTool
label: Uses new pan organ Azimuth tool to annotate cell types for any organ

hints:
  DockerRequirement:
    dockerPull: sennet/rna-data-products-azimuth
baseCommand: /opt/pan_organ_azimuth.py


inputs:
  raw_h5ad_file:
    type: File
    inputBinding:
      position: 0
  organism:
    type: string?
    inputBinding:
      position: 1
  tissue:
    type: string?
    inputBinding:
      position: 2


outputs:
  annotated_raw_h5ad_file:
    type: File
    outputBinding:
      glob: '*_raw.*'
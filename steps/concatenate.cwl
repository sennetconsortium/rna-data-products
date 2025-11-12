cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: sennet/rna-data-products-python
baseCommand: /opt/annotate_concatenate.py

inputs:
  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?
    inputBinding:
      position: 0

  data_directory:
    type: Directory
    doc: Paths to directory containing processed datasets
    inputBinding:
      position: 2

  uuids_file:
    type: File
    doc: Path to file containing uuids of all datasets to be indexed
    inputBinding:
      position: 3

  tissue:
    type: string?
    doc: string description of tissue
    inputBinding:
      position: 4

outputs:
  raw_h5ad_file:
    type: File
    outputBinding:
      glob: "*_raw.h5ad"
    doc: h5ad file containing raw expression data

  data_product_metadata:
    type: File
    outputBinding: 
      glob: "*.json"
    doc: metadata for web app
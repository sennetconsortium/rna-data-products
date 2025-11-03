cwlVersion: v1.0
class: CommandLineTool
label: Perform secondary analysis on raw data product

hints:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-sketching
baseCommand: /opt/sketching.py

inputs: 
    partially_processed_h5ad_file:
        type: File
        doc: The raw h5ad file
        inputBinding:
            position: 0

    tissue:
        type: string?
        doc: optional tissue type
        inputBinding:
            position: 1

outputs:
    sketched_h5ad_file:
        type: File
        outputBinding:
            glob: "*_processed.h5ad"
        doc: h5ad file with some secondary analysis processing

cwlVersion: v1.0
class: CommandLineTool
label: Perform secondary analysis on raw data product

hints:
  DockerRequirement:
    dockerPull: sennet/rna-data-products-python
baseCommand: /opt/secondary_analysis_pt1.py

inputs: 
    annotated_raw_h5ad_file:
        type: File
        doc: The raw h5ad file
        inputBinding:
            position: 0

    tissue:
        type: string?
        doc: optional tissue type
        inputBinding:
            position: 1
    
    uuids_file:
        type: File
        doc: File with UUIDs and patient metadata
        inputBinding:
            position: 2
    
    data_product_metadata:
        type: File
        doc: data product metadata
        inputBinding: 
            position: 3

outputs:
    final_raw_h5mu_file:
        type: File
        outputBinding:
            glob: "*_raw.h5mu"
        doc: annotated h5mu file with additional obs columns

    partially_processed_h5ad_file:
        type: File
        outputBinding:
            glob: "*_processed.h5ad"
        doc: h5ad file with secondary analysis processing
    
    updated_data_product_metadata:
        type: File
        outputBinding:
            glob: "*.json"
        doc: data product metadata

cwlVersion: v1.0
class: CommandLineTool
label: Perform secondary analysis on raw data product

hints:
  DockerRequirement:
    dockerPull: sennet/rna-data-products-python
baseCommand: /opt/secondary_analysis_pt2.py

inputs: 
    sketched_h5ad_file:
        type: File
        doc: The raw h5ad file
        inputBinding:
            position: 0

    tissue:
        type: string?
        doc: optional tissue type
        inputBinding:
            position: 1
    
    updated_data_product_metadata:
        type: File
        doc: data product metadata
        inputBinding: 
            position: 2

outputs:
    processed_h5ad_file:
        type: File
        outputBinding:
            glob: "*_processed.h5ad"
        doc: h5ad file with secondary analysis processing
    
    processed_h5mu_file:
        type: File
        outputBinding:
            glob: "*_processed.h5mu"
        doc: mudata object for EPIC registration

    umap_png:
        type: File
        outputBinding:
            glob: "*.png"
        doc: umap png
    
    final_data_product_metadata:
        type: File
        outputBinding:
            glob: "*.json"
        doc: final data product metadata with all cell type counts

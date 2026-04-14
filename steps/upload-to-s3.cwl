cwlVersion: v1.0
class: CommandLineTool
label: Uploads the annotated and concatenated h5ad files and the umap png

hints:
  DockerRequirement:
    dockerPull: sennet/rna-data-products-python
baseCommand: /opt/upload_to_s3.py

inputs:
    final_raw_h5mu_file:
        type: File
        doc: The raw h5mu file
        inputBinding:
            position: 0

    processed_h5mu_file:
        type: File
        doc: The processed h5mu file
        inputBinding:
            position: 1

    umap_png:
        type: File
        doc: PNG of UMAP
        inputBinding:
            position: 2

    final_data_product_metadata:
        type: File
        doc: data product metadata json
        inputBinding:
            position: 3

    shinycell_dir:
        type: Directory
        doc: directory containing shinyapp files
        inputBinding:
            position: 4

    tissue:
        type: string
        doc: tissue type
        inputBinding:
            position: 5

    access_key_id:
        type: string
        doc: AWS access key id
        inputBinding:
            position: 6

    secret_access_key:
        type: string
        doc: AWS secret access key
        inputBinding:
            position: 7

outputs:
    finished_text:
        type: File
        outputBinding:
            glob: "*.txt"
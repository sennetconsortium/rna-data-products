#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

requirements:
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:

  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?

  data_directory:
    label: "Path to directory containing processed RNA datasets"
    type: Directory

  uuids_file:
    label: "Path to a file containing a list of uuids for the dataset to be indexed"
    type: File

  tissue:
    label: "String description of tissue type"
    type: string?
  
  organism:
    label: "human or mouse, default is human"
    type: string?
    default: "human"

  access_key_id:
    label: "AWS access key id"
    type: string

  secret_access_key:
    label: "AWS secret access key"
    type: string

outputs:

  final_raw_h5mu_file:
    outputSource: secondary-analysis/final_raw_h5mu_file
    type: File
  processed_h5mu_file:
    outputSource: secondary-analysis/processed_h5mu_file
    type: File
  umap_png:
    outputSource: secondary-analysis/umap_png
    type: File
  final_data_product_metadata:
    outputSource: secondary-analysis/final_data_product_metadata
    type: File
  shinycell_dir:
    outputSource: make-shinycell/shinycell_dir
    type: Directory
  deepscence_continuous_plot:
    outputSource: secondary-analysis/deepscence_continuous_plot
    type: File
    label: "Umap with coloring by DeepScence scores"
  deepscence_binary_plot:
    outputSource: secondary-analysis/deepscence_binary_plot
    type: File
    label: "Umap with coloring by DeepScence binary results"
  deepscence_plot:
    outputSource: deepscence/deepscence_plot
    type: File
    
steps:

  - id: concatenate
    in:
      - id: enable_manhole
        source: enable_manhole
      - id: data_directory
        source: data_directory
      - id: uuids_file
        source: uuids_file
      - id: tissue
        source: tissue

    out:
      - raw_h5ad_file
      - data_product_metadata
    run: steps/concatenate.cwl
    label: "Concatenates h5ad data files in directory"

  - id: deepscence
    in:
      - id: raw_h5ad_file
        source: concatenate/raw_h5ad_file
    out:
      - raw_h5ad_with_ds
      - deepscence_plot
    run: steps/deepscence.cwl

  - id: secondary-analysis
    in:
      - id: raw_h5ad_with_ds
        source: deepscence/raw_h5ad_with_ds
      - id: tissue
        source: tissue
      - id: uuids_file
        source: uuids_file
      - id: data_product_metadata
        source: concatenate/data_product_metadata

    out:
      - final_raw_h5mu_file
      - processed_h5ad_file
      - processed_h5mu_file
      - umap_png
      - final_data_product_metadata
      - deepscence_continuous_plot
      - deepscence_binary_plot
    run: steps/secondary-analysis.cwl
    label: "Runs secondary anaylsis on annotated and concatenated data and annotates with pan organ Azimuth"

  - id: make-shinycell
    in:
      - id: processed_h5ad_file
        source: secondary-analysis/processed_h5ad_file
      - id: tissue
        source: tissue
      - id: metadata_file
        source: secondary-analysis/final_data_product_metadata
    out:
      - shinycell_dir

    run: steps/make_shinycell.cwl
    label: "Creates the shiny cell app for the data product"

  - id: upload-to-s3
    in:
      - id: final_raw_h5mu_file
        source: secondary-analysis/final_raw_h5mu_file
      - id: processed_h5mu_file
        source: secondary-analysis/processed_h5mu_file
      - id: umap_png
        source: secondary-analysis/umap_png
      - id: final_data_product_metadata
        source: secondary-analysis/final_data_product_metadata
      - id: shinycell_dir
        source: make-shinycell/shinycell_dir
      - id: tissue
        source: tissue
      - id: access_key_id
        source: access_key_id
      - id: secret_access_key
        source: secret_access_key

    out:
      - finished_text
    run: steps/upload-to-s3.cwl
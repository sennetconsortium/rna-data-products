#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

requirements:
  ScatterFeatureRequirement: {}

inputs:

  raw_h5ad_file:
    label: "Concatenated h5ad file"
    type: File

  tissue:
    label: "Two letter tissue code"
    type: string?

  uuids_file:
    label: "Path to a file containing a list of uuids for the dataset to be indexed"
    type: File

  data_product_metadata:
    label: "Data product metadata JSON file"
    type: File

outputs:

  final_raw_h5mu_file:
    outputSource: secondary-analysis-pt1/final_raw_h5mu_file
    type: File
  processed_h5ad_file:
    outputSource: secondary-analysis-pt2/processed_h5ad_file
    type: File
  processed_h5mu_file:
    outputSource: secondary-analysis-pt2/processed_h5mu_file
    type: File
  umap_png:
    outputSource: secondary-analysis-pt2/umap_png
    type: File
  final_data_product_metadata:
    outputSource: secondary-analysis-pt2/final_data_product_metadata
    type: File

steps:

  - id: secondary-analysis-pt1
    in:
      - id: raw_h5ad_file
        source: raw_h5ad_file
      - id: tissue
        source: tissue
      - id: uuids_file
        source: uuids_file
      - id: data_product_metadata
        source: data_product_metadata

    out:
      - final_raw_h5mu_file
      - partially_processed_h5ad_file
      - updated_data_product_metadata
    run: secondary-analysis/secondary-analysis-pt1.cwl
    label: "Runs secondary anaylsis on concatenated data"

  - id: azimuth-annotate
    in:
      - id: partially_processed_h5ad_file
        source: secondary-analysis-pt1/partially_processed_h5ad_file
      - id: tissue
        source: tissue
    out:
      - annotated_h5ad_file
    run: secondary-analysis/azimuth-annotate.cwl

  - id: sketching
    in:
      - id: annotated_h5ad_file
        source: azimuth-annotate/annotated_h5ad_file
    out:
      - sketched_h5ad_file
    run: secondary-analysis/sketching.cwl
    label: sketches cells based on PCA density

  - id: secondary-analysis-pt2
    in:
      - id: sketched_h5ad_file
        source: sketching/sketched_h5ad_file
      - id: tissue
        source: tissue
      - id: updated_data_product_metadata
        source: secondary-analysis-pt1/updated_data_product_metadata

    out:
      - processed_h5ad_file
      - processed_h5mu_file
      - umap_png
      - final_data_product_metadata
    run: secondary-analysis/secondary-analysis-pt2.cwl
    label: "Runs secondary anaylsis on annotated and concatenated data"

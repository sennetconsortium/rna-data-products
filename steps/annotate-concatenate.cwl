#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

requirements:
  ScatterFeatureRequirement: {}

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
    type: string?
    default: "human"

outputs:

  annotated_raw_h5ad_file:
    outputSource: azimuth-annotate/annotated_raw_h5ad_file
    type: File
  data_product_metadata:
    outputSource: concatenate/data_product_metadata
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
    run: annotate-concatenate/concatenate.cwl
    label: "Concatenates h5ad data files in directory"

  - id: azimuth-annotate
    in: 
      - id: raw_h5ad_file
        source: concatenate/raw_h5ad_file
      - id: tissue
        source: tissue
      - id: organism
        source: organism
    
    out:
      - annotated_raw_h5ad_file
    run: annotate-concatenate/azimuth-annotate.cwl
    label: "Runs azimuth on the file created in the previous step"
    
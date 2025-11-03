#!/usr/bin/env python3

from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import anndata
import json
import mudata as md
import numpy as np
import os
import pandas as pd
import scanpy as sc

def add_patient_metadata(obs, uuids_df):
    merged = uuids_df.merge(obs, left_on="uuid", right_on="dataset", how="inner")
    merged = merged.set_index(obs.index)
    merged = merged.drop(columns=["Unnamed: 0"], errors='ignore')
    merged = merged.fillna(np.nan)
    merged["age"] = pd.to_numeric(merged["age"], errors='coerce')
    obs = obs.loc[:, ~obs.columns.str.contains("^Unnamed")]
    return merged


def add_file_sizes(data_product_metadata, raw_size):
    data_product_metadata["Raw File Size"] = raw_size


def add_cell_counts(data_product_metadata, cell_counts, total_cell_count):
    with open(data_product_metadata, "r") as json_file:
        metadata = json.load(json_file)
    metadata["Raw Cell Type Counts"] = cell_counts
    metadata["Raw Total Cell Count"] = total_cell_count
    return metadata


def main(
    raw_h5ad_file: Path,
    uuids_tsv: Path,
    data_product_metadata: Path,
    tissue: str = None,
):
    raw_output_file_name = f"{tissue}_raw.h5mu" if tissue else "rna_raw.h5mu"
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "rna_processed.h5ad"
    )

    adata = anndata.read_h5ad(raw_h5ad_file)
    dataset_info = pd.read_csv(uuids_tsv, sep="\t")
    annotated_obs = add_patient_metadata(adata.obs, dataset_info)
    total_cell_count = adata.obs.shape[0]
    cell_type_counts = adata.obs["final_level_labels"].value_counts().to_dict()
    adata.uns["cell_type_counts"] = json.dumps(cell_type_counts)
    metadata = add_cell_counts(
        data_product_metadata, cell_type_counts, total_cell_count
    )
    adata.obs = annotated_obs

    with open(data_product_metadata, "r") as infile:
        metadata = json.load(infile)
    uuid = metadata["Data Product UUID"]
    # Convert to MuData and add Obj x Analyte requirements
    adata.obs['object_type'] = 'cell'
    adata.uns['analye_class'] = 'RNA'
    adata.uns['protocol'] = 'https://github.com/hubmapconsortium/rna-data-products.git'
    if 'azimuth_label' in adata.obs_keys():
        azimuth = adata.obs[['full_hierarchical_labels', 'final_level_labels', 'final_level_confidence', 'full_consistent_hierarchy', 'azimuth_broad', 'azimuth_medium', 'azimuth_fine', 'CL_Label', 'CL_ID']]
        adata.obsm['annotation'] = pd.DataFrame(adata.obs['final_level_labels'])
        adata.obsm['azimuth'] = azimuth
        adata.uns['azimuth'] = {
            'label': 'Cell Ontology Annotation',
            'ontologyID': 'predicted_CLID',
            'mechanism': 'machine',
            'protocol': "10.1016/j.cell.2021.04.048",
        }
    mdata = md.MuData({f"{uuid}_raw": adata})
    mdata.uns['epic_type '] = ['analyses', 'annotations']
    mdata.write(raw_output_file_name)

    raw_file_size = os.path.getsize(raw_output_file_name)
    add_file_sizes(metadata, raw_file_size)
    with open(f"{uuid}.json", "w") as outfile:
        json.dump(metadata, outfile)

    print("Processing data product")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["unscaled"] = adata.X.copy()
    sc.pp.scale(adata, max_value=10)

    sc.pp.pca(adata, n_comps=50)
    sc.tl.embedding_density(adata, basis="pca")  

    adata.write(processed_output_file_name)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("uuids_file")
    p.add_argument("data_product_metadata")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole
        manhole.install(activate_on="USR1")

    main(args.raw_h5ad_file, args.uuids_file, args.data_product_metadata, args.tissue)

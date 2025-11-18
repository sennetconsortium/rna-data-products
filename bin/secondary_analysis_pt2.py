#!/usr/bin/env python3

from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import anndata
import json
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import os
import pandas as pd
import scanpy as sc
import uuid


def add_cell_counts(data_product_metadata, cell_counts, total_cell_count):
    with open(data_product_metadata, "r") as json_file:
        metadata = json.load(json_file)
    metadata["Processed Cell Type Counts"] = cell_counts
    metadata["Processed Total Cell Count"] = total_cell_count
    return metadata


def add_file_sizes(data_product_metadata, processed_size):
    data_product_metadata["Processed File Size"] = processed_size
    uuid = data_product_metadata["Data Product UUID"]
    with open(f"{uuid}.json", "w") as outfile:
        json.dump(data_product_metadata, outfile)


def main(h5ad_file: Path, data_product_metadata: Path, tissue: str=None):
    adata = anndata.read_h5ad(h5ad_file)
    processed_output_file_name = (
        f"{tissue}_processed" if tissue else "rna_processed"
    )
    total_cell_count = adata.obs.shape[0]
    with open(data_product_metadata, "r") as infile:
        metadata = json.load(infile)
    uuid = metadata["Data Product UUID"]
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)
    sc.tl.umap(adata)

    # leiden clustering
    sc.tl.leiden(adata)
    sc.tl.rank_genes_groups(adata, "leiden")

    if "CL_Label" in adata.obs_keys():

        non_na_values = adata.obs.final_level_labels.dropna()
        counts_dict = non_na_values.value_counts().to_dict()
        keep_cell_types = [
            cell_type for cell_type in counts_dict if counts_dict[cell_type] > 1
        ]
        adata_filter = adata[adata.obs.final_level_labels.isin(keep_cell_types)]
        sc.tl.rank_genes_groups(
            adata_filter, "CL_Label", key_added="rank_genes_groups_cell_types"
        )
        adata.uns = adata_filter.uns

    if "CL_Label" in adata.obs_keys():
        cell_type_counts = adata.obs["CL_Label"].value_counts().to_dict()
        adata.uns["cell_type_counts"] = json.dumps(cell_type_counts)
        metadata = add_cell_counts(
            data_product_metadata, cell_type_counts, total_cell_count
        )
    else:
        cell_type_counts = {}
        metadata = add_cell_counts(
            data_product_metadata, cell_type_counts, total_cell_count
        )

    with plt.rc_context():
        sc.pl.umap(adata, color="leiden", show=False)
        plt.savefig(f"{uuid}.png")

    # Convert to MuData and add Obj x Analyte requirements
    if 'annotation' in adata.obsm_keys():
        adata.obsm['annotation']['leiden'] = adata.obs['leiden']
    else:
        adata.obsm['annotation'] = pd.DataFrame(adata.obs['leiden'])
    adata.obsm['leiden'] = pd.DataFrame(adata.obs['leiden'])
    adata.uns['leiden'] = {
        'label': 'Leiden Clusters',
        'mechanism': 'machine',
        'protocol': "10.1186/s13059-017-1382-0",
    }
    if 'CL_Label' in adata.obs_keys():
        azimuth = adata.obs[['full_hierarchical_labels', 'final_level_labels', 'final_level_confidence', 'full_consistent_hierarchy', 'azimuth_broad', 'azimuth_medium', 'azimuth_fine', 'CL_Label', 'CL_ID']]
        adata.obsm['annotation'] = pd.DataFrame(adata.obs['CL_Label'])
        adata.obsm['azimuth'] = azimuth
        adata.uns['azimuth'] = {
            'label': 'Cell Ontology Annotation',
            'ontologyID': 'predicted_CLID',
            'mechanism': 'machine',
            'protocol': "10.1016/j.cell.2021.04.048",
        }
    mdata = mu.MuData({f"{uuid}_processed": adata})
    mdata.uns["epic_type "] = ['analyses', 'annotations']

    print(f"Writing {processed_output_file_name}")
    adata.write(f"{processed_output_file_name}.h5ad")
    mdata.write(f"{processed_output_file_name}.h5mu")
    processed_file_size = os.path.getsize(f"{processed_output_file_name}.h5mu")
    add_file_sizes(metadata, processed_file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("data_product_metadata", type=Path)
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.raw_h5ad_file, args.data_product_metadata, args.tissue)
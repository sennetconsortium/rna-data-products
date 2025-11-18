#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, check_call
from typing import Optional
import anndata
import pandas as pd
import scanpy as sc
import squidpy as sq
import muon as mu

import panhumanpy
from plot_utils import new_plot
import matplotlib.pyplot as plt
import json

CLID_MAPPING = "/opt/pan-human-azimuth-crosswalk.csv"
model_version = "v1"


def map_to_clid(adata_obs: pd.DataFrame):
    reference = pd.read_csv(CLID_MAPPING, header=10)
    label_to_cl_label = dict(zip(reference['Annotation_Label'], reference['CL_Label']))
    label_to_cl_id = dict(zip(reference['Annotation_Label'], reference['CL_ID']))

    adata_obs['CL_Label'] = adata_obs['final_level_labels'].map(label_to_cl_label)
    adata_obs['CL_ID'] = adata_obs['final_level_labels'].map(label_to_cl_id)

    return adata_obs


def main(
        secondary_analysis_matrix: Path,
        tissue: str = None,
):
    adata = anndata.read_h5ad(secondary_analysis_matrix)
    outfile = f'{tissue}_processed.h5ad' if tissue else 'rna_processed.h5ad'

    if "unscaled" in adata.layers:
        adata.X = adata.layers["unscaled"]

    adata = adata[:, ~adata.var.hugo_symbol.isna()]

    adata.write('secondary_analysis_hugo.h5ad')
    azimuth_annotate_command = f"annotate secondary_analysis_hugo.h5ad -fn hugo_symbol -mv {model_version}"
    check_call(azimuth_annotate_command, shell=True)
    ct_adata = anndata.read_h5ad('secondary_analysis_hugo_ANN.h5ad')
    secondary_analysis_adata = anndata.AnnData(X=adata.X, var=adata.var, obs=ct_adata.obs,
                                               obsm=ct_adata.obsm, uns=ct_adata.uns)

    for key in adata.uns.keys():
        secondary_analysis_adata.uns[key] = adata.uns[key]

    secondary_analysis_adata.obs = map_to_clid(secondary_analysis_adata.obs)
    secondary_analysis_adata.uns["pan_human_azimuth_crosswalk"] = {
        "title": "Cell type annotations for pan-human Azimuth, v1.0",
        "description": (
            "This crosswalk maps cell type annotations from pan-human Azimuth to the "
            "Cell Ontology (Version IRI: http://purl.obolibrary.org/obo/cl/releases/2025-04-10/cl.owl)."
        ),
        "url": "https://cdn.humanatlas.io/digital-objects/ctann/pan-human-azimuth/latest/assets/pan-human-azimuth-crosswalk.csv",
        "publisher": "HuBMAP",
        "creators": ["Aleix Puig-Barbe"],
        "project_lead": "Katy Börner",
        "reviewers": ["Bruce W. Herr II", "Katy Börner"],
        "processor": "HRA Digital Object Processor, v0.9.0",
        "date_published": "2025-06-15",
        "date_last_processed": "2025-06-12",
        "funders": [
            "National Institutes of Health (OT2OD033756)",
            "National Institutes of Health (OT2OD026671)"
        ],
        "license": "CC BY 4.0",
        "dashboard": "https://apps.humanatlas.io/dashboard/data"
    }
    secondary_analysis_adata.uns["annotation_metadata"] = {}
    secondary_analysis_adata.uns["annotation_metadata"]["is_annotated"] = True
    secondary_analysis_adata.uns["annotation_metadata"]["model_version"] = model_version
    secondary_analysis_adata.uns["annotation_metadata"]["panhumanpy_version"] = panhumanpy.__version__

    for key in adata.obsm:
        secondary_analysis_adata.obsm[key] = adata.obsm[key]

    for key in adata.layers:
        secondary_analysis_adata.layers[key] = adata.layers[key]

    secondary_analysis_adata.write(outfile)

    calculated_metadata_dict = {"annotation_tools": ["Pan-human Azimuth"], "object_types": ["CL:0000000"]}
    with open('calculated_metadata.json', 'w') as f:
        json.dump(calculated_metadata_dict, f)

    cell_type_manifest_dict = {}

    for column_header in ['azimuth_broad', 'azimuth_medium', 'azimuth_fine', 'final_level_labels', 'CL_ID']:
        sub_dict = {
            val: int((secondary_analysis_adata.obs[column_header] == val).sum())
            for val in secondary_analysis_adata.obs[column_header].unique()
        }
        # Remove NaN key if it exists
        sub_dict = {k: v for k, v in sub_dict.items() if not pd.isna(k)}
        cell_type_manifest_dict[column_header] = sub_dict

    with open('cell_type_manifest.json', 'w') as f:
        json.dump(cell_type_manifest_dict, f)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('raw_h5ad_file', type=Path)
    p.add_argument('tissue', type=str, nargs="?")
    args = p.parse_args()

    main(
        args.raw_h5ad_file,
        args.tissue,
    )

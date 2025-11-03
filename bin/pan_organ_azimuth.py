#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from subprocess import run, check_call
from typing import Optional
import anndata
import pandas as pd
import scanpy as sc
import squidpy as sq
import mudata as md
from plot_utils import new_plot
import matplotlib.pyplot as plt
import json

CLID_MAPPING = "/opt/pan-human-azimuth-crosswalk.csv"


def map_to_clid(adata_obs: pd.DataFrame):
    
    reference = pd.read_csv(CLID_MAPPING, header=10)
    label_to_cl_label = dict(zip(reference['Annotation_Label'], reference['CL_Label']))
    label_to_cl_id = dict(zip(reference['Annotation_Label'], reference['CL_ID']))

    adata_obs['CL_Label'] = adata_obs['final_level_labels'].map(label_to_cl_label)
    adata_obs['CL_ID'] = adata_obs['final_level_labels'].map(label_to_cl_id)

    return adata_obs


def main(
        raw_h5ad_file: Path,
        organism: str,
        tissue: str = None,
):
    if organism == "human":
        if raw_h5ad_file.suffix == ".h5mu":
            mudata = md.read_h5mu(raw_h5ad_file)
            adata = mudata.mod["rna"]
        else:
            adata = anndata.read_h5ad(raw_h5ad_file)

        if "unscaled" in adata.layers:
            adata.X = adata.layers["unscaled"]

        adata = adata[:,~adata.var.hugo_symbol.isna()]

        adata.write('raw_hugo.h5ad')
        azimuth_annotate_command = f"annotate raw_hugo.h5ad -fn hugo_symbol"
        check_call(azimuth_annotate_command, shell=True)
        ct_adata = anndata.read_h5ad('raw_hugo_ANN.h5ad')
        annotated_adata = anndata.AnnData(X=adata.X, var=adata.var, obs=ct_adata.obs,
                                                   obsm = ct_adata.obsm, uns=ct_adata.uns)
        print(annotated_adata)

        for key in adata.uns.keys():
            annotated_adata.uns[key] = adata.uns[key]

        annotated_adata.obs = map_to_clid(annotated_adata.obs)
        annotated_adata.uns["pan_human_azimuth_crosswalk"] = {
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
        annotated_adata.uns["annotation_metadata"] = {}
        annotated_adata.uns["annotation_metadata"]["is_annotated"] = True

        for key in adata.obsm:
            annotated_adata.obsm[key] = adata.obsm[key]

        for key in adata.layers:
            annotated_adata.layers[key] = adata.layers[key]

        if raw_h5ad_file.suffix == ".h5mu":
            mudata.mod["rna"] = annotated_adata
            mudata.write_h5mu(f"{tissue}_raw.h5mu" if tissue else "RNA_raw.h5mu")
        else:
            annotated_adata.write(f"{tissue}_raw.h5ad" if tissue else "RNA_raw.h5ad")

        calculated_metadata_dict = {"annotation_tools": ["Pan-human Azimuth"], "object_types": ["CL:0000000"]}
        with open('calculated_metadata.json', 'w') as f:
            json.dump(calculated_metadata_dict, f)

        cell_type_manifest_dict = {}

        for column_header in ['azimuth_broad', 'azimuth_medium', 'azimuth_fine', 'final_level_labels', 'CL_ID']:
            sub_dict = {
                val: int((annotated_adata.obs[column_header] == val).sum())
                for val in annotated_adata.obs[column_header].unique()
            }
            # Remove NaN key if it exists
            sub_dict = {k: v for k, v in sub_dict.items() if not pd.isna(k)}
            cell_type_manifest_dict[column_header] = sub_dict

        with open('cell_type_manifest.json', 'w') as f:
            json.dump(cell_type_manifest_dict, f)

    else:
        if raw_h5ad_file.suffix == ".h5mu":
            mudata = md.read_h5mu(raw_h5ad_file)
            mudata.uns["annotation_metadata"] = {}
            mudata.uns["annotation_metadata"]["is_annotated"] = False
            mudata.write_h5mu(f"{tissue}_raw.h5mu" if tissue else "RNA_raw.h5mu")
        else:
            adata = anndata.read_h5ad(raw_h5ad_file)
            adata.uns["annotation_metadata"] = {}
            adata.uns["annotation_metadata"]["is_annotated"] = False
            adata.write(f"{tissue}_raw.h5ad" if tissue else "RNA_raw.h5ad")

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('raw_h5ad_file', type=Path)
    p.add_argument('organism', type=str)
    p.add_argument('tissue', type=str, nargs="?")
    args = p.parse_args()

    main(
        args.raw_h5ad_file,
        args.organism,
        args.tissue,
    )
    
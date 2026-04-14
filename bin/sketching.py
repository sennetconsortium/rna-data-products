#!/usr/bin/env python3

from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import anndata
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import sclkme
import uuid


def main(h5ad_file: Path, tissue:str=None):
    # keep cells with density of less than 0.95 to start, possibly go down to 0.9
    adata = anndata.read_h5ad(h5ad_file)
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "rna_processed.h5ad"
    )
    low_density_cells = adata.obs["pca_density"] < 0.95
    num_cells = low_density_cells.sum()
    sclkme.tl.sketch(adata, n_sketch=num_cells, use_rep="X_pca", method="gs", key_added="gs")
    adata_sketch_gs = adata[adata.obs["gs_sketch"]]
    adata_sketch_gs.write(processed_output_file_name)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.raw_h5ad_file, args.tissue)
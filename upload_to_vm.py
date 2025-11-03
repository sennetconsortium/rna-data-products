#!/usr/bin/env python3
import json
import os
from argparse import ArgumentParser
from pathlib import Path


def get_uuid(data_product_metadata):
    with open(data_product_metadata, "r") as json_file:
        metadata = json.load(json_file)
    uuid = metadata["Data Product UUID"]
    return uuid


def upload_to_vm(umap_png, metadata_json, shiny_cell_dir, uuid):
    os.system(
        f"scp {umap_png} /opt/repositories/vm024-dev/data-products-ui/pipeline_outputs/{uuid}.png"
    )
    os.system(
        f"scp {metadata_json} /opt/repositories/vm024-dev/data-products-ui/pipeline_outputs/{uuid}.json"
    )
    os.system(
        f"mkdir /opt/repositories/vm024-dev/data-products-ui/shiny-server/{uuid}"
    )
    os.system(
        f"scp {shiny_cell_dir}/* /opt/repositories/vm024-dev/data-products-ui/shiny-server/{uuid}"
    )


def main(umap_png, metadata_json, shiny_cell_dir):
    uuid = get_uuid(metadata_json)
    upload_to_vm(umap_png, metadata_json, shiny_cell_dir, uuid)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("umap_png", type=Path)
    p.add_argument("data_product_metadata", type=Path)
    p.add_argument("shiny_cell_dir", type=Path)
    args = p.parse_args()

    main(args.umap_png, args.data_product_metadata, args.shiny_cell_dir)

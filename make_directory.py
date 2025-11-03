#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import requests
import subprocess


def main(uuids_file: Path, tissue: str):
    uuids = pd.read_csv(uuids_file, sep="\t")["uuid"].dropna()
    h5ads_base_directory = Path(f"{tissue}_h5ads")
    h5ads_base_directory.mkdir(exist_ok=True)

    for uuid in uuids:
        h5ads_directory = h5ads_base_directory / uuid
        h5ads_directory.mkdir(parents=True, exist_ok=True)
        out_file = h5ads_directory / "expr.h5ad"
        url = f"https://assets.api.sennetconsortium.org/{uuid}/expr.h5ad"

        # Check if file exists before downloading
        try:
            response = requests.head(url, allow_redirects=True)
            if response.status_code != 200:
                print(f"File for {uuid} does not exist on server (status {response.status_code}).")
                continue
        except requests.RequestException as e:
            print(f"Error checking {uuid}: {e}")
            continue
        print(f"Downloading {uuid} ...")
        try:
            subprocess.run(["wget", "-q", url, "-O", str(out_file)], check=True)
        except subprocess.CalledProcessError:
            print(f"Failed to download {uuid}")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("uuids_file", type=Path)
    p.add_argument("tissue", type=str)

    args = p.parse_args()

    main(args.uuids_file, args.tissue)

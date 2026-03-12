#!/usr/bin/env python3
from argparse import ArgumentParser

import pandas as pd
import requests
import yaml
import json

organ_types_yaml_file = "bin/organ_types.yaml"
organ_uberon_file = "bin/organs.json"


def get_uuids(organ_uberon: str, organism: str):
    params = {
        "status": "Published",
        "dataset_type": "RNAseq [Salmon]",
    }

    if organ_uberon:
        params["origin_samples.organ"] = organ_uberon

    url = "https://search.api.sennetconsortium.org/param-search/datasets"
    response = requests.get(url, params=params)
    print("Request URL:", response.url)
    print("Response status code: ", response.status_code)

    # Handle a successful response
    if response.status_code == 200:
        return process_response(response, organism)

    # Handle 303 redirection
    elif response.status_code == 303:
        redirection_url = (
            response.text.strip()
        )  # Get redirection URL from response body
        print("Following redirection URL: ", redirection_url)

        # Make a request to the redirection URL
        redirection_response = requests.get(redirection_url)
        if redirection_response.status_code == 200:
            return process_response(redirection_response, organism)

    # Handle other error responses
    else:
        print(f"Error {response.status_code}: {response.text}")
        return [], [], []


def process_response(response, organism):
    data = response.json()
    items = data

    uuids = []
    sennet_ids = []
    donor_metadata_list = []

    for item in items:
        # skip bulk datasets
        dataset_info = item.get("dataset_info", "")
        if dataset_info and "bulk" in dataset_info.lower():
            continue
        # extract info from 'sources'
        sources = item.get("sources", [])
        if not sources:
            continue
        source = sources[0]
        if source.get("source_type", "").lower() == organism.lower():
            uuids.append(item.get("uuid"))
            sennet_ids.append(item.get("sennet_id"))
        # Attempt to extract donor metadata, if available
            metadata = source.get("metadata")
            if organism == "mouse":
                print(metadata)
                donor_metadata_list.append(extract_mouse_metadata(metadata))
            elif metadata.get("living_donor_data"):
                donor_metadata = metadata.get("living_donor_data")
                donor_metadata_list.append(extract_donor_metadata(donor_metadata))
            else:
                donor_metadata = metadata.get("organ_donor_data")
                donor_metadata_list.append(extract_donor_metadata(donor_metadata))

    return (
        uuids,
        sennet_ids,
        donor_metadata_list,
    )


def extract_mouse_metadata(donor_metadata):
    mouse_info = {
        "bedding": donor_metadata.get("bedding"),
        "cage_enhancements": donor_metadata.get("cage_enhancements"),
        "data_of_birth_or_fertilization": donor_metadata.get("data_of_birth_or_fertilization"),
        "date_of_death": donor_metadata.get("date_of_death"),
        "diet": donor_metadata.get("diet"),
        "euthanization_method": donor_metadata.get("euthanization_method"),
        "is_deceased": donor_metadata.get("is_deceased"),
        "is_embryo": donor_metadata.get("is_embryo"),
        "light_cycle": donor_metadata.get("light_cycle"),
        "local_lifespan_data": donor_metadata.get("local_lifespan_data"),
        "rack_setup": donor_metadata.get("rack_setup"),
        "room_health_status": donor_metadata.get("room_health_status"),
        "room_temperature": donor_metadata.get("room_temperature"),
        "sex": donor_metadata.get("sex"),
        "strain": donor_metadata.get("strain"),
        "strain_rrid": donor_metadata.get("strain_rrid"),
        "water_source": donor_metadata.get("water_source"),
    }
    return mouse_info


def extract_donor_metadata(donor_metadata):
    donor_info = {
        "age": None,
        "sex": None,
        "height": None,
        "weight": None,
        "bmi": None,
        "cause_of_death": None,
        "race": None,
        "medical_history": None,
        "abo_blood_type": None,
        "mechanism_of_injury": None,
    }
    for item in donor_metadata:
        if item.get("grouping_concept_preferred_term") == "ABO blood group system":
            donor_info["abo_blood_type"] = item.get("data_value")
        elif item.get("grouping_concept_preferred_term") == "Age":
            donor_info["age"] = item.get("data_value") + " " + item.get("units")
        elif item.get("grouping_concept_preferred_term") == "Body Mass Index":
            donor_info["bmi"] = item.get("data_value") + " " + item.get("units")
        elif item.get("grouping_concept_preferred_term") == "Cause of Death":
            donor_info["cause_of_death"] = item.get("data_value")
        elif item.get("grouping_concept_preferred_term") == "Height":
            donor_info["height"] = item.get("data_value") + " " + item.get("units")
        elif item.get("grouping_concept_preferred_term") == "Mechanism of Injury":
            donor_info["mechanism_of_injury"] = item.get("data_value")
        elif item.get("grouping_concept_preferred_term") == "Race":
            donor_info["race"] = item.get("data_value")
        elif item.get("grouping_concept_preferred_term") == "Sex":
            donor_info["sex"] = item.get("data_value")
        elif item.get("grouping_concept_preferred_term") == "Medical History":
            donor_info["medical_history"] = item.get("data_value")
        elif item.get("grouping_concept_preferred_term") == "Weight":
            donor_info["weight"] = item.get("data_value") + " " + item.get("units")
    return donor_info


def get_organ_uberon(organ_name):
    term = organ_name.lower().strip()
    with open(organ_uberon_file) as f:
        data = json.load(f)
    for entry in data:
        if entry.get("term", "").lower() == term:
            return entry["organ_uberon"]
    for entry in data:
        cat = entry.get("category")
        if cat and cat.get("term", "").lower() == term:
            return cat["organ_uberon"]
    return None


def main(tissue_type: str, organism:str):
    organ_dict = yaml.load(open(organ_types_yaml_file), Loader=yaml.BaseLoader)
    for key in organ_dict:
        organ_dict[key] = organ_dict[key]["description"]
    uberon_code = get_organ_uberon(tissue_type)
    uuids_list, sennet_ids_list, donor_metadata = get_uuids(uberon_code, organism)
    uuids_df = pd.DataFrame()
    uuids_df["uuid"] = pd.Series(uuids_list, dtype=str)
    uuids_df["sennet_id"] = pd.Series(sennet_ids_list, dtype=str)
    donor_metadata_df = pd.DataFrame(donor_metadata)
    result_df = pd.concat([uuids_df, donor_metadata_df], axis=1)
    key_for_tissue = [key for key, value in organ_dict.items() if value == tissue_type]
    if key_for_tissue:
        output_file_name = f"{key_for_tissue[0].lower()}_{organism}.tsv"
    else:
        output_file_name = "rna.tsv"
    result_df['organism'] = organism
    result_df['tissue'] = tissue_type
    print(result_df)
    result_df.to_csv(output_file_name, sep="\t")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("tissue_type", type=str, nargs="?", help="Type of tissue (optional)")
    p.add_argument("organism", type=str)
    args = p.parse_args()

    main(args.tissue_type, args.organism)

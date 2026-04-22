# rna-data-products
A Python and CWL pipeline for concatenating HuBMAP RNA-seq [Salmon] data into data products per organ and one large RNA-seq [Salmon] data product.
## Pipeline steps
* Create a UUIDs TSV file with all UUIDs and HuBMAP IDs of public processed data wanted for the run.
* With the UUIDs TSV, create a data directory of all H5ADs needed for the run.
* Make an AWS access key id and a secret access key to upload the files to S3 bucket.
* Annotate and concatenate a raw data product and a processed data product.
* Upload the UMAP and data product metadata to VM
## Requirements
Check the list of python packages in `docker/requirements.txt`
## How to run
### Step 1
`python3 make_uuids_tsv.py [tissue_type]`
### Step 2
`python3 make_directory.py /hive/hubmap/data/ [uuids_file] [tissue_type]`
### Step 3 
`cwltool pipeline.cwl --[data_directory] --[uuids_file] --[tissue_type] --[access_key_id] --[secret_access_key]`
### Step 4
`python3 upload_to_ec2.py [umap_png] [data_product_metadata] [shiny_cell_dir] [ssh_key]`

# Rosetta data

We provide raw Rosetta data as well as processed Rosetta datasets that have duplicates, outliers, and NaN values removed. The data is hosted on Zenodo. 

## Raw Rosetta data
Raw Rosetta data comes in the form of SQLite databases. 
There are separate databases for each of the local datasets as well as the global dataset.

| Dataset         | Filename           | Compressed size (approx.) | Uncompressed size (approx.) | Direct download | MD5 checksum                       |
|-----------------|--------------------|---------------------------|-----------------------------|-----------------|------------------------------------|
| GFP             | avgfp.tar.gz       | 5 GB                      | 15 GB                       | Link            | `039141a2693c5e3907ff34cf19df8ee0` |
| DLG4            | dlg4.tar.gz        | 6 GB                      | 15 GB                       | Link            | `3c734d96bc636a477338e3b638f44e9f` |
| GB1             | gb1.tar.gz         | 3 GB                      | 8 GB                        | Link            | `ca4d2e90e81fd2a0d5f00017e5f0b145` |
| GB1-IgG Binding | gb1_binding.tar.gz | 1 GB                      | 4 GB                        | Link            | `59443943105b662adfbdb77dceb04a1e` |
| GRB2            | grb2.tar.gz        | 5 GB                      | 16 GB                       | Link            | `f0c30f43bf4ce47deb0423ffd00fefc8` |
| Pab1            | pab1.tar.gz        | 5 GB                      | 13 GB                       | Link            | `284e7bf4351814cb72a92c5c18a43de1` |
| TEM-1           | tem-1.tar.gz       | 5 GB                      | 15 GB                       | Link            | `e0295dd42447b2795771963afe201867` |
| Ube4b           | ube4b.tar.gz       | 5 GB                      | 13 GB                       | Link            | `bbefcd91863e18c1ac72af8f1c5f7b06` |
| Global          | global.tar.gz      | 9 GB                      | 21 GB                       | Link            | `74aad528e294433b0c55e0b8c5b75213` |

Note the GB1-IgG binding raw data only contains the binding scores, whereas the processed GB1-IgG binding dataset listed below contains both the binding and standard scores. 
The processed dataset was created by combining the raw GB1-IgG binding data with the raw GB1 standard data.

The script [parse_rosetta_data.py](https://github.com/gitter-lab/metl/blob/main/code/parse_rosetta_data.py) and notebook [generate_rosetta_dataset.ipynb](https://github.com/gitter-lab/metl/blob/main/notebooks/generate_rosetta_dataset.ipynb), which are located in the [metl](https://github.com/gitter-lab/metl) repository, show how to process these raw databases into the processed datasets listed below.

If needed, you can load this data directly into a dataframe with the Pandas [pd.read_sql()](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_sql.html) function.
However, we found the built-in Pandas function to be slow with large databases. Instead, we recommend the [ConnectorX](https://github.com/sfu-db/connector-x) package to speed up loading.

```python
import connectorx as cx

db_fn = "path/to/database.db"
conn = f"sqlite://{db_fn}"

query = "SELECT * FROM `variant`"
result = cx.read_sql(conn, query)
```

Additionally, the HDF5 format is substantially faster to load and can be used with Pandas. If you need to regularly load these raw databases, it may be worthwhile to convert them to HDF5. 

## Processed Rosetta datasets
Each processed Rosetta dataset has its own directory containing the following:
- The dataset in three formats (.tsv, .db, and .h5 files), all containing the same data
- A list of PDB files corresponding to the variants in the dataset (pdb_fns.txt)
- A splits directory containing train, validation, and test splits we used for pretraining
- Standardization parameters computed on the train set (in the splits directory)

Processed Rosetta datasets can be used directly with the main [metl](https://github.com/gitter-lab/metl) repository to pretrain models.

The HDF5 file uses the Pandas "fixed" format and can be loaded into a dataframe using:
```python
pd.read_hdf(fn, key="variant")
```

| Dataset         | Filename           | Compressed size (approx.) | Uncompressed size (approx.) | Direct download | MD5 checksum                       |
|-----------------|--------------------|---------------------------|-----------------------------|-----------------|------------------------------------|
| GFP             | avgfp.tar.gz       | 12 GB                     | 34 GB                       | Link            | `a5536c91289cca054ad4de07fe8494e0` |
| DLG4            | dlg4-2022.tar.gz   | 14 GB                     | 40 GB                       | Link            | `81677115c1318e7720ebf3b866443a81` |
| GB1             | gb1.tar.gz         | 7 GB                      | 21 GB                       | Link            | `2c6efa9bc2d6a8f3e8b801f62397907a` |
| GB1-IgG Binding | gb1_binding.tar.gz | 4 GB                      | 13 GB                       | Link            | `68760620b893be9e8aacb15c7df4f01b` |
| GRB2            | grb2.tar.gz        | 12 GB                     | 38 GB                       | Link            | `6eeb625c4c934873808cf4d2a11fcf3e` |
| Pab1            | pab1.tar.gz        | 12 GB                     | 33 GB                       | Link            | `aa384c9984a9d8499b8468f2108e8c0e` |
| TEM-1           | tem-1.tar.gz       | 13 GB                     | 37 GB                       | Link            | `78a91362359a06070917529b41b82a63` |
| Ube4b           | ube4b.tar.gz       | 12 GB                     | 34 GB                       | Link            | `be071466be3c08fe7fec2431ea404b91` |
| Global          | global.tar.gz      | 20 GB                     | 51 GB                       | Link            | `02221f43363c23b06b156900fbd86957` |

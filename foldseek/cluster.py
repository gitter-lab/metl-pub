import pandas as pd
import yaml
from Bio import SeqIO

def parse_datasets_yaml():
    # Load YAML file
    file_path = 'data/dms_data/datasets.yml'
    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)
    return data

def make_clusters(df):

    final_dict = {}
    for row in df.itertuples():
        final_dict[row.cluster]=[]

    for row in df.itertuples():
        final_dict[row.cluster].append(row.pdb)

    return final_dict


def return_multiple_pdb_clusters(cluster_dict):
    multiple_pdb_cluster ={}
    for key in cluster_dict.keys():
        if len(cluster_dict[key]) > 1:
            multiple_pdb_cluster[key]=cluster_dict[key]

    return multiple_pdb_cluster

def check_if_clusters_overlap_with_test_dms(multiple_pdb_cluster,datasets,fasta_dict):

    for dataset_key in datasets.keys():
        dataset_pdb_id = datasets[dataset_key]['rosettafy_pdb_fn'].split('.')[0]
        for pdb_key in multiple_pdb_cluster.keys():
            if dataset_pdb_id in multiple_pdb_cluster[pdb_key]:
                print(f'\nfound cluster for dms id: {dataset_key}')
                print(f'>{dataset_pdb_id}\n{datasets[dataset_key]["wt_aa"]}')
                for cluster_pdb_id in multiple_pdb_cluster[pdb_key]:
                    if cluster_pdb_id==dataset_pdb_id:
                        pass
                    else:
                        print(f">{cluster_pdb_id}\n{fasta_dict[cluster_pdb_id]}")


def parse_fasta(fasta_fn):

    # Load sequences from the FASTA file
    seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_fn, "fasta")}

    return seq_dict
def main(fn,fasta_fn):
    print('\n')
    print(100*'-')
    print(f'cluster:{fn}')
    datasets= parse_datasets_yaml()
    fasta_dict= parse_fasta(fasta_fn)
    df= pd.read_csv(fn,sep='\t')

    new_columns = ['cluster','pdb']
    # Save current headers as a row
    current_headers = pd.DataFrame([df.columns], columns=new_columns)

    # Rename columns
    df.columns = new_columns

    # Concatenate the new row with the DataFrame
    df = pd.concat([current_headers, df], ignore_index=True)

    all_clusters =  make_clusters(df)

    multiple_clusters_only=return_multiple_pdb_clusters(all_clusters)



    check_if_clusters_overlap_with_test_dms(multiple_clusters_only, datasets,fasta_dict)




if __name__ == '__main__':
    for file in ['foldseek/metl_0.9_cluster.tsv','foldseek/metl_mmseqs_0.9_cluster.tsv']:
        main(file,'foldseek/metl_all_seqs.fasta')



### SCRIPT for filtering out the pairs from the AHoJ-DB batch output

import os
import pandas as pd
import math
import sys
sys.path.append('..')
import pymol_utils
sys.path.append('../B-create-dataset')
import filter_utils
import zipfile
import json

# only holo structures with good resolution (max 4 Angstroms) and without any apo structures
# OUTPUT_PATH = '/home/skrhakv/CryptoBench/data/A-filter-ahojdb/holo_only_label_seq_mapping'

# all holo structures, those with apo structures as well
# OUTPUT_PATH = '/home/skrhakv/CryptoBench/data/A-filter-ahojdb/holo_only_with_apo'

# all holos with filtered ligands (filter_utils.filter_valid_ligands(...))
# OUTPUT_PATH = '/home/skrhakv/CryptoBench/data/A-filter-ahojdb/all_holo'

# this is the same but it doesn't include other ligands from the HOLO structure (query_pocket_info.csv)
# OUTPUT_PATH = '/home/skrhakv/CryptoBench/data/A-filter-ahojdb/old_all_apoholo_all_ligands'
 
OUTPUT_PATH = '/home/skrhakv/CryptoBench/data/A-filter-ahojdb/all_apoholo_all_ligands'
INPUT_PATH = '/home/skrhakv/CryptoBench/data/ahoj-db/ahojdb_v2c/data'
LIGANDS_PATH = '/home/skrhakv/CryptoBench/data/B-create-dataset/holo_only_data'


def get_label_seq_binding_site(holo_info, batch, job):
    """Get the binding site in label_seq_id format from the residue_mappings.zip zipped file."""
    # unzip residue mapping
    zip_file_path = f'{INPUT_PATH}/{batch}/{job}/residue_mappings.zip'
    extract_dir = f'{INPUT_PATH}/{batch}/{job}/residue_mappings'
    residue_mapping_file = f'{extract_dir}/combined_simplified_pocket_info.csv'

    if not os.path.exists(residue_mapping_file):
        if os.path.exists(zip_file_path):
            with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
        else:
            return None

    # load residue mapping    
    if not os.path.exists(residue_mapping_file):
        return None

    pdb_id = holo_info["structure"].iloc[0]
    chain = holo_info["chains3"].iloc[0]
    residue_mapping = pd.read_csv(residue_mapping_file)
    residue_mapping = residue_mapping[(residue_mapping["query_struct"] == pdb_id) &
                                        (residue_mapping["POI"] == '_'.join(job.split('-')[1:]))]['PDBe_index']

    # store the residues in label_seq_id format
    mapped_residues = f'{pdb_id}_{chain} and ( resi {"+".join(map(str, residue_mapping))} )'         

    return mapped_residues

def main():
    output_file_location = f'{OUTPUT_PATH}/pairs.csv'

    for batch in os.listdir(INPUT_PATH):
        for job in os.listdir(f'{INPUT_PATH}/{batch}'):
            
            # skip log files
            if job == 'log':
                continue
            
            print(f'processing {job}, batch: {batch}')

            if not os.path.exists(f'{INPUT_PATH}/{batch}/{job}/query_pocket_info.csv') or not os.path.exists(f'{INPUT_PATH}/{batch}/{job}/pocket_selections.csv'):
                continue
            # this happens in the data and I don't know why, rather skip it

            if not os.path.exists(f'{INPUT_PATH}/{batch}/{job}/apo_filtered_sorted_results.csv'):
                continue

            # skip if apo structures exist
            apo_info = pd.read_csv(
                f'{INPUT_PATH}/{batch}/{job}/apo_filtered_sorted_results.csv')

            # load holo info
            holo_info = pd.read_csv(
                f'{INPUT_PATH}/{batch}/{job}/query_pocket_info.csv')
            apo_info['resolution'] = apo_info['resolution'].replace('-', math.inf)

            # consider holo when there is complete pocket and good resolution
            holo_info['resolution'] = holo_info['resolution'].replace(
            '-', math.inf)

            # this shouldn't happen, but just in case
            if holo_info.empty:
                continue

            # skip if holo structure is multi-chain
            if '-' in str(holo_info['chains3'].iloc[0]):
                continue 

            assert holo_info.shape[0] == 1, f'Error: more than one holo structure found for {job}, batch: {batch}'

            apo_headers = apo_info.to_dict()
            holo_headers = holo_info.to_dict()
            apo_headers = {f'apo_{k}': [] for k, v in apo_headers.items()}
            holo_headers = {f'holo_{k}': [] for k, v in holo_headers.items()}
            if 'apo_RoG_distance' not in apo_headers:
                apo_headers['apo_RoG_distance'] = []

            apo_holo_pairs = apo_headers | holo_headers
            apo_holo_pairs['apo_pocket_selection'] = []
            apo_holo_pairs['holo_pocket_selection'] = []

            pocket_selections = pd.read_csv(
                f'{INPUT_PATH}/{batch}/{job}/pocket_selections.csv', header=None)

            holo_pocket_selection = pymol_utils.get_pocket_selection(
                pocket_selections, holo_info.iloc[0])

            if apo_info.empty:
                apo_holo_pairs['apo_pocket_selection'].append(None)
                apo_holo_pairs['holo_pocket_selection'].append(
                    holo_pocket_selection)
                for header, column in holo_info.iloc[0].items():
                    apo_holo_pairs[f'holo_{header}'].append(column)
                for header in apo_headers.keys():
                    if header != 'apo_pocket_selection':
                        apo_holo_pairs[header].append(None)
                
            else:
                for apo_row1 in apo_info.itertuples():
                    apo_pocket_selection1 = pymol_utils.get_pocket_selection(
                        pocket_selections, apo_row1)
                    if apo_pocket_selection1 == None:
                        continue
                    apo_holo_pairs[f'apo_pocket_selection'].append(
                        apo_pocket_selection1)
                    apo_holo_pairs[f'holo_pocket_selection'].append(
                        holo_pocket_selection)
                    for header, column in apo_info.loc[apo_row1.Index].items():
                        apo_holo_pairs[f'apo_{header}'].append(column)
                    for header, column in holo_info.iloc[0].items():
                        apo_holo_pairs[f'holo_{header}'].append(column)

            df_apo_holo_pairs = pd.DataFrame.from_dict(apo_holo_pairs)

            df_apo_holo_pairs.to_csv(output_file_location, mode='a',
                                 header=not os.path.exists(output_file_location))

if __name__ == '__main__':
    main()

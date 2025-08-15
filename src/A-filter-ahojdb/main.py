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
# with commented 66-68th line
OUTPUT_PATH = '/home/skrhakv/CryptoBench/data/A-filter-ahojdb/all_holo'
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

            # uncommented in data/A-filter-ahojdb/holo_only_with_apo:
            #
            # this happens in the data and I don't know why, rather skip it
            # if not os.path.exists(f'{INPUT_PATH}/{batch}/{job}/apo_filtered_sorted_results.csv'):
            #     continue
            if not os.path.exists(f'{INPUT_PATH}/{batch}/{job}/query_pocket_info.csv') or not os.path.exists(f'{INPUT_PATH}/{batch}/{job}/pocket_selections.csv'):
                continue

            # load holo info
            holo_info = pd.read_csv(
                f'{INPUT_PATH}/{batch}/{job}/query_pocket_info.csv')

            # consider holo when there is complete pocket and good resolution
            holo_info['resolution'] = holo_info['resolution'].replace(
            '-', math.inf)

            # skip if no holo structures with good resolution exist
            if holo_info.empty:
                continue
            # skip if holo structure is multi-chain
            if '-' in str(holo_info['chains3'].iloc[0]):
                continue 

            assert holo_info.shape[0] == 1, f'Error: more than one holo structure found for {job}, batch: {batch}'

            # retrieve appropriate pocket selections (PDB and AF) from the pocket_selections.csv file
            pocket_selections = pd.read_csv(
                f'{INPUT_PATH}/{batch}/{job}/pocket_selections.csv', header=None)
            pdb_pocket_selection = pymol_utils.get_pocket_selection(
                pocket_selections, holo_info.iloc[0])
            alphafold_pocket_selection = pymol_utils.get_alphafold_pocket_selection(
                pocket_selections, holo_info.iloc[0])
            if alphafold_pocket_selection is None:
                print(f'Error: no AF pocket selection found for {job}, batch: {batch}')
                continue
            # add the selection to the holo_info dataframe
            holo_info['pdb_pocket_selection'] = [pdb_pocket_selection]
            holo_info['alphafold_pocket_selection'] = [alphafold_pocket_selection]

            mapped_binding_site = get_label_seq_binding_site(holo_info, batch, job)
            
            if mapped_binding_site is None:
                continue
            holo_info['label_seq_pdb_pocket_selection'] = [mapped_binding_site]

            # skip if ligand is not valid
            holo_info = filter_utils.filter_valid_ligands(
                holo_info, LIGANDS_PATH, query_poi_name='query_POI')
            if holo_info.empty:
                continue

            # append it to the output file
            holo_info.to_csv(output_file_location, mode='a',
                header=not os.path.exists(output_file_location))
    
if __name__ == '__main__':
    main()

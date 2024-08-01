import pandas as pd
import sys
import json
import pickle
import shutil

sys.path.append('../B-create-dataset')
import filter_utils
import pocket
import pocket_selection_parser

import filter_utils

LIGANDS_DATA_PATH = '../../data/B-create-dataset/ahoj-v2'
DATASET_PATH = '../../data/C-remove-holo-homomers/ahoj-v2'
AHOJ_DB_PATH = '../../data/A-filter-ahojdb-v2/pairs.csv'
OUTPUT_PATH = '/home/vit/Projects/cryptobench/data/E-add-noncryptic-pockets/ahoj-v2'
NON_CRYPTIC_DS_PATH = f'{OUTPUT_PATH}/cryptobench/auxiliary-data/non-cryptic-pockets'

def find_additional_noncryptic_pockets(df, apo_pdb_id, apo_chain_id):
    additional_noncryptic = df[(df['apo_structure'] == apo_pdb_id) & (df['apo_chains3'] == apo_chain_id)]
    # print('get valid pairs ...')
    additional_noncryptic = filter_utils.get_well_defined_pairs(additional_noncryptic)
    additional_noncryptic = additional_noncryptic[(additional_noncryptic['apo_pocket_rms'] < filter_utils.MAX_POCKET_RMSD)]
    
    # print('filter non-valid ligands out ...')
    if len(additional_noncryptic) < 1:
        return None
    return filter_utils.filter_valid_ligands(additional_noncryptic, LIGANDS_DATA_PATH)

def concat_pd(pd1, pd2):
    if pd1 is None and pd2 is None:
        return None
    elif pd1 is not None and pd2 is None:
        return pd1
    elif pd2 is None and pd1 is not None:
        return pd2
    else:
        return pd.concat([pd1, pd2])


def check_final_dataset():
    """Sanity check that every key in noncryptic-pockets.json is present in dataset.json"""
    with open(f'{NON_CRYPTIC_DS_PATH}/noncryptic-pockets.json') as f:
        noncryptic = json.load(f)

    with open(f'{DATASET_PATH}/dataset.json') as f:
        ds = json.load(f)

    for key in noncryptic.keys():
        assert key in ds

def main():
    with open(f'{DATASET_PATH}/dataset.json') as f:
        ds = json.load(f)

    # get noncryptic pockets
    additional_pockets = {}
    # file is too big - load in chunks
    for chunk_iteratior, chunk in enumerate(pd.read_csv(AHOJ_DB_PATH, chunksize=2000000)):
        for apo_iterator, (apo_pdb_id, holo_structures) in enumerate(ds.items()):  
            apo_chain = holo_structures[0]['apo_chain']
            apo_id = apo_pdb_id + apo_chain
            
            print(f'Processing {apo_id}; chunk {chunk_iteratior}, {apo_iterator} / {len(ds.keys())}')

            if apo_id not in additional_pockets:
                additional_pockets[apo_id] = None
            
            # deal with potential multichain structures - search for each chain alone
            for apo_subchain in apo_chain.split('-'):
                new_pockets = find_additional_noncryptic_pockets(chunk, apo_pdb_id, apo_subchain)
                additional_pockets[apo_id] = concat_pd(additional_pockets[apo_id], new_pockets)

            # search additional pockets for the whole multichain
            if '-' in apo_chain:
                new_pockets = find_additional_noncryptic_pockets(chunk, apo_pdb_id, apo_chain)
                additional_pockets[apo_id] = concat_pd(additional_pockets[apo_id], new_pockets)       
    
    # the previous step takes a long time (~40min) -> cache the intermediate results just in case
    with open(f'{OUTPUT_PATH}/additional_pockets.pickle', 'wb') as handle:
        pickle.dump(additional_pockets, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # parse the extracted pockets - use existing code from filter_utils 
    parser = pocket_selection_parser.pocket_selection_parser()
    parsed_additional_pockets = {}
    for apo_pdb_id, holo_structures in additional_pockets.items():
        id = (apo_pdb_id[:4], apo_pdb_id[4:])
        if holo_structures is not None and len(holo_structures) > 0:
            # idk why there were NaN values in holo_pocket_selection -> filter them out
            holo_structures = holo_structures[holo_structures['holo_pocket_selection'].notna()]
            _, pockets = parser.parse_pocket_selections(holo_structures)
            parsed_additional_pockets[id] = pockets
        else: parsed_additional_pockets[id] = []

    # save the dataset
    filter_utils.save_dataset(parsed_additional_pockets, NON_CRYPTIC_DS_PATH, 'noncryptic-pockets.json')
    check_final_dataset()

    # copy README.md to the cryptobench folder
    shutil.copy('../README.md', f'{OUTPUT_PATH}/cryptobench')
    
if __name__ == '__main__':
    main()
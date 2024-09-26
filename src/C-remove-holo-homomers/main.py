import json
import itertools
import biotite.database.rcsb as rcsb
import os
from biotite.structure.io.pdbx import get_sequence
import biotite.structure.io.pdbx as pdbx

CIF_FILES_PATH = '/home/vit/Projects/deeplife-project/data/cif_files'
COMPUTE_FOR_CRYPTIC = False

if COMPUTE_FOR_CRYPTIC:
    OUTPUT_PATH = '/home/vit/Projects/cryptobench/data/C-remove-holo-homomers/ahoj-v2'
    INPUT_PATH = '/home/vit/Projects/cryptobench/data/B-create-dataset/ahoj-v2'
else:
    OUTPUT_PATH = '/home/vit/Projects/cryptobench/data/E-add-noncryptic-pockets/ahoj-v2/cryptobench/auxiliary-data/non-cryptic-pockets'
    INPUT_PATH = OUTPUT_PATH

def do_correction(dataset):
    """Some Ids got outdated. Let's updated them"""
    for apo_pdb_id, holo_structures in dataset.items():
        for i, holo_structure in enumerate(holo_structures):
            if holo_structure['uniprot_id'] == 'O36607':
                dataset[apo_pdb_id][i]['uniprot_id'] = 'Q2HRB6'
    return dataset

def main():
    if COMPUTE_FOR_CRYPTIC:
        file = 'dataset.json'
    else:
        file = 'noncryptic-pockets.json'
    
    with open(f'{INPUT_PATH}/{file}') as f:
        dataset = json.load(f)

    dataset = do_correction(dataset)

    dataset_without_homomers = {}
    counter = 0
    identified_homomers = 0
    for apo_pdb_id, holo_structures in dataset.items():
        dataset_without_homomers[apo_pdb_id] = []
        for _, holo_group in itertools.groupby(holo_structures, key=lambda row: [set(row['apo_pocket_selection']), row['holo_pdb_id'], row['ligand']]):
            holo_group = list(holo_group)

            if len(holo_group) > 1:

                # pdb id is the same for the whole group
                holo_pdb_id = holo_group[0]['holo_pdb_id']

                print(holo_pdb_id)

                # get chains from the group:
                group_chains = []
                for holo_structure in holo_group:
                    group_chains.append(holo_structure['holo_chain'].split('-'))

                # load sequences
                cif_file_path = rcsb.fetch(
                    holo_pdb_id, "cif", CIF_FILES_PATH)
                cif_file = pdbx.CIFFile.read(cif_file_path)
                structure_sequences = [i.as_array() for i in list(
                    cif_file[holo_pdb_id.upper()]['entity_poly'].values())]
                structure_sequences_chains = structure_sequences[6]

                # add the first item because we are sure that it doesn't have homomeric duplicate as there is no such structure yet
                dataset_without_homomers[apo_pdb_id].append(holo_group[0])
                # if any pair of group_chains items has the same sequence then it is homomeric and one of them can be deleted
                is_homomeric_idx = []
                for i, pair_item1 in enumerate(group_chains):
                    for ii, pair_item2 in enumerate(group_chains):
                        # these pairs were already checked, skip those:
                        if ii <= i:
                            continue
                        # these pairs were already identified as homomeric, also skip:
                        if i in is_homomeric_idx or ii in is_homomeric_idx:
                            continue

                        # check each position in multichain structures
                        is_every_chain_in_multichain_homomeric = True
                        for multichain_index in range(len(pair_item1)):
                            is_this_position_homomeric = False
                            # check each homomer set whether each chain at particular position belongs to that homomer
                            for homomer_chains in structure_sequences_chains: 
                                homomer_chains = set(homomer_chains.split(','))
                                if pair_item1[multichain_index] in homomer_chains and pair_item2[multichain_index] in homomer_chains:
                                    is_this_position_homomeric = True

                            # the position was not homomeric
                            if not is_this_position_homomeric:
                                is_every_chain_in_multichain_homomeric = False
                                break
                        
                        # this record doesn't have (yet) homomeric duplicate
                        if not is_every_chain_in_multichain_homomeric:            
                            dataset_without_homomers[apo_pdb_id].append(holo_group[ii])
                        else: 
                            is_homomeric_idx.append(ii)

                print(f'Group chains: {group_chains}')
                print(f'Homomers: {structure_sequences_chains}')
                print(f'Identified homomers: {[group_chains[i] for i in is_homomeric_idx]}')
                identified_homomers += len([group_chains[i] for i in is_homomeric_idx])
            else:
                dataset_without_homomers[apo_pdb_id].append(holo_group[0])

        counter += 1

        if counter % 100:
            with open(f'{OUTPUT_PATH}/{file}', 'w', encoding='utf-8') as f:
                json.dump(dataset_without_homomers, f, ensure_ascii=False, indent=4)        
    print(f'total of identified homomers: {identified_homomers}')

if __name__ == '__main__':
    main()

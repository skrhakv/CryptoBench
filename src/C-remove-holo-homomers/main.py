import json
import itertools
import biotite.database.rcsb as rcsb
import os
from biotite.structure.io.pdbx import get_sequence
import biotite.structure.io.pdbx as pdbx

CIF_FILES_PATH = '/home/vit/Projects/deeplife-project/data/cif_files'
OUTPUT_PATH = '/home/vit/Projects/cryptobench/data/C-remove-holo-homomers/ahoj-v2'
INPUT_PATH = '/home/vit/Projects/cryptobench/data/B-create-dataset/ahoj-v2'


def main():
    with open(f'{INPUT_PATH}/dataset.json') as f:
        dataset = json.load(f)

    dataset_without_homomers = {}
    counter = 0
    for apo_pdb_id, holo_structures in dataset.items():
        holo_ids_to_be_removed = set()
        for _, holo_group in itertools.groupby(holo_structures, key=lambda row: [set(row['apo_pocket_selection']), row['holo_pdb_id'], row['ligand']]):
            holo_group = list(holo_group)

            if len(holo_group) > 1:

                # pdb id is the same for the whole group
                holo_pdb_id = holo_group[0]['holo_pdb_id']

                print(holo_pdb_id)

                # get chains from the group:
                group_chains = []
                for holo_structure in holo_group:
                    group_chains.append(holo_structure['holo_chain'])

                # load sequences
                cif_file_path = rcsb.fetch(
                    holo_pdb_id, "cif", CIF_FILES_PATH)
                cif_file = pdbx.CIFFile.read(cif_file_path)
                structure_sequences = [i.as_array() for i in list(
                    cif_file[holo_pdb_id.upper()]['entity_poly'].values())]
                structure_sequences_chains = structure_sequences[6]

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

                        for homomer_chains in structure_sequences_chains: 
                            homomer_chains = set(homomer_chains.split(','))
                            if pair_item1 in homomer_chains and pair_item2 in homomer_chains:
                                is_homomeric_idx.append(ii)

                print(f'Group chains: {group_chains}')
                print(f'Homomers: {structure_sequences_chains}')
                print(f'Identified homomers: {[group_chains[i] for i in is_homomeric_idx]}')

                for i in is_homomeric_idx:
                    holo_ids_to_be_removed.add(f'{holo_pdb_id}{group_chains[i]}')
        dataset_without_homomers[apo_pdb_id] = []
        for holo_structure in dataset[apo_pdb_id]:
            if f'{holo_structure["holo_pdb_id"]}{holo_structure["holo_chain"]}' in holo_ids_to_be_removed:
                continue
            else:
                dataset_without_homomers[apo_pdb_id].append(holo_structure)

        counter += 1

        if counter % 100:
            with open(f'{OUTPUT_PATH}/dataset.json', 'w', encoding='utf-8') as f:
                json.dump(dataset_without_homomers, f, ensure_ascii=False, indent=4)

        


if __name__ == '__main__':
    main()

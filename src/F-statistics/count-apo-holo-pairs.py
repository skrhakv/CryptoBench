#!/usr/bin/env python3

import json
import pandas as pd
import matplotlib.pyplot as plt

INPUT_PATH = '/home/vit/Projects/cryptobench/data/C-remove-holo-homomers/rigid-dataset'
ALL_PAIRS_CSV_PATH = '/home/vit/Projects/cryptobench/data/B-create-dataset/rigid-dataset/filtered_rmsd.csv'


with open(f'{INPUT_PATH}/dataset.json') as f:
    dataset = json.load(f)

all_pairs = pd.read_csv(ALL_PAIRS_CSV_PATH)

# the number of alternative apo-holo pairs for each apo structure - min(len(apo_pairs), len(holo_pairs))
number_of_pairs = {}

for apo_pdb_id in dataset.keys():
    holo_structure = dataset[apo_pdb_id][0]
    apo_query_POI = f"{holo_structure['ligand_chain']}_{holo_structure['ligand']}_{holo_structure['ligand_index']}"

    out = all_pairs[(all_pairs['apo_query_POI'] == apo_query_POI) & (all_pairs['apo_structure'] == apo_pdb_id)].iloc[0]
    minimum = min(out['number_of_alternative_holos'], out['number_of_alternative_apos'])
    # minimum = out['number_of_alternative_holos'] + out['number_of_alternative_apos']
    if minimum not in number_of_pairs:
        number_of_pairs[minimum] = 0
    number_of_pairs[minimum] += 1

for k in sorted(number_of_pairs.keys()):
    print(f'{k}: {number_of_pairs[k]}')

plt.bar(list(number_of_pairs.keys()), number_of_pairs.values(), color='g')
plt.show()

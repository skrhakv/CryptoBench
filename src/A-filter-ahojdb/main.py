### SCRIPT for filtering out the pairs from the AHoJ-DB batch output

import os
import pandas as pd
import math
import sys
sys.path.append('..')
import pymol_utils
sys.path.append('../B-create-dataset')
import filter_utils

MINIMAL_RESOLUTION = 2.5
OUTPUT_PATH = '/home/skrhakv/CryptoBench/data/A-filter-ahojdb/holo_only_data'
INPUT_PATH = '/home/skrhakv/CryptoBench/data/ahoj-db/ahojdb_v2c/data'
LIGANDS_PATH = '/home/skrhakv/CryptoBench/data/B-create-dataset/holo_only_data'

def main():
    output_file_location = f'{OUTPUT_PATH}/pairs.csv'

    for batch in os.listdir(INPUT_PATH):
        for job in os.listdir(f'{INPUT_PATH}/{batch}'):
            
            # skip log files
            if job == 'log':
                continue
            
            print(f'processing {job}, batch: {batch}')

            # this happens in the data and I don't know why, rather skip it
            if not os.path.exists(f'{INPUT_PATH}/{batch}/{job}/apo_filtered_sorted_results.csv'):
                continue

            # skip if apo structures exist
            apo_info = pd.read_csv(
                f'{INPUT_PATH}/{batch}/{job}/apo_filtered_sorted_results.csv')
            
            if not apo_info.empty:
                continue

            # load holo info
            holo_info = pd.read_csv(
                f'{INPUT_PATH}/{batch}/{job}/query_pocket_info.csv')

            # consider holo when there is complete pocket and good resolution
            holo_info['resolution'] = holo_info['resolution'].replace(
            '-', math.inf)
            holo_info = holo_info[(holo_info['resolution'] < MINIMAL_RESOLUTION)]

            # skip if no holo structures with good resolution exist
            if holo_info.empty:
                continue
            # skip if holo structure is multi-chain
            if '-' in str(holo_info['chains3'].iloc[0]):
                continue 

            assert holo_info.shape[0] == 1, f'Error: more than one holo structure found for {job}, batch: {batch}'

            # retrieve appropriate pocket selection from the pocket_selections.csv file
            pocket_selections = pd.read_csv(
                f'{INPUT_PATH}/{batch}/{job}/pocket_selections.csv', header=None)
            pocket_selection = pymol_utils.get_pocket_selection(
                pocket_selections, holo_info.iloc[0])

            # add the selection to the holo_info dataframe
            holo_info['pocket_selection'] = [pocket_selection]

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

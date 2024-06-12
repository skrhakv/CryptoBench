### SCRIPT for filtering out the pairs from the AHoJ-DB batch output

import os
import time
import pandas as pd
import math
import sys
sys.path.append('..')
import pymol_utils
OUTPUT_PATH = 'output/pairs.csv'
INPUT_PATH = '/home/vit/Projects/ahoj2-extraction/src/get-pairs/tmp-data/'

MINIMAL_RESOLUTION = 2.5

def main(batch_id, input_path, output_path):
    output_file_location = f'{output_path}/pairs.csv'

    start = time.time()
    counter = 1
    for filename in os.listdir(f'{input_path}/{batch_id}-pairs'):
        print(
            f'processing {filename}: {counter}/{len(os.listdir(f"{input_path}/{batch_id}-pairs"))}, batch: {batch_id}')
        counter += 1

        os.makedirs(f'output', exist_ok=True)
        holo_info = pd.read_csv(
            f'{input_path}/{batch_id}-pairs/{filename}/query_pocket_info.csv')
        pocket_selections = pd.read_csv(
            f'{input_path}/{batch_id}-pairs/{filename}/pocket_selections.csv', header=None)

        headers = holo_info.to_dict()
        d1 = {f'apo_{k}': [] for k, v in headers.items()}
        d2 = {f'holo_{k}': [] for k, v in headers.items()}
        apo_holo_pairs = d1 | d2
        apo_holo_pairs['apo_pocket_selection'] = []
        apo_holo_pairs['holo_pocket_selection'] = []

        # consider holo when there is complete pocket and good resolution
        holo_info['resolution'] = holo_info['resolution'].replace(
            '-', math.inf)
        holo_info = holo_info[(holo_info['resolution'] < MINIMAL_RESOLUTION)]

        if holo_info.empty:
            continue
        
        pocket_size = holo_info['#resis3'].iloc[0]

        apo_info = pd.read_csv(
            f'{input_path}/{batch_id}-pairs/{filename}/apo_filtered_sorted_results.csv')
        apo_info['resolution'] = apo_info['resolution'].replace('-', math.inf)

        # consider apo where there is complete pocket and good resolution
        apo_info = apo_info[(apo_info['resolution'] < MINIMAL_RESOLUTION) &
                            (apo_info['pocket_len'] == pocket_size)]
        
        holo_pocket_selection = pymol_utils.get_pocket_selection(
            pocket_selections, holo_info.iloc[0])
        
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
    end = time.time()
    print(end - start)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])

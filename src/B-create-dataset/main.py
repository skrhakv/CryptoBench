import pandas as pd
import filter_utils
import os

SINGLE_CHAIN_ONLY = False

INPUT_PATH = '../../data/A-filter-ahojdb'
INPUT_CSV = f'{INPUT_PATH}/pairs.csv'
OUTPUT_PATH = '../../data/B-create-dataset'
if SINGLE_CHAIN_ONLY:
    OUTPUT_PATH += '/single-chain-with-covalent-ligands'
else:
    OUTPUT_PATH += '/multi-chain-with-covalent-ligands'
PREFILTERED_CSV = f'{OUTPUT_PATH}/filtered_rmsd.csv'
DATASET_DATAFRAME = f'{OUTPUT_PATH}/dataset_dataframe.csv'

# DEFINE FILTER THRESHOLDS
MAX_POCKET_RMSD = 2
MIN_TM_SCORE = 0.5
MAX_POCKET_DISTANCE = 4
MIN_ROG_PERCENTAGE = 0.8
MAX_ROG_PERCENTAGE = 1.2


def load_ahojdb():
    print('Load AHoJ-DB ...')

    if not os.path.exists(PREFILTERED_CSV):
        ahojdb_df = pd.read_csv(INPUT_CSV)
        filtered_rmsd_df = ahojdb_df[(
            ahojdb_df['apo_pocket_rms'] > MAX_POCKET_RMSD)]
        filtered_rmsd_df.to_csv(PREFILTERED_CSV)
        del ahojdb_df
    else:
        filtered_rmsd_df = pd.read_csv(PREFILTERED_CSV)

    if SINGLE_CHAIN_ONLY:
        filtered_rmsd_df = filter_utils.remove_multichain_pockets(
            filtered_rmsd_df)
        
    # filter AHoJ-DB according to the defined thresholds
    filtered_rmsd_df = \
        filtered_rmsd_df[(filtered_rmsd_df['apo_pocket_rms'] > MAX_POCKET_RMSD) &
                         (filtered_rmsd_df['apo_tm_score'] > MIN_TM_SCORE) &
                         (filtered_rmsd_df['apo_tm_score_i'] > MIN_TM_SCORE) &
                         (filtered_rmsd_df['apo_pocket_dist'] < MAX_POCKET_DISTANCE) &
                         # (filtered_rmsd_df['holo_nonquery_ligs'].isna()) &
                         (filtered_rmsd_df['apo_RoG'] * MAX_ROG_PERCENTAGE > filtered_rmsd_df['holo_RoG']) &
                         (filtered_rmsd_df['apo_RoG'] * MIN_ROG_PERCENTAGE < filtered_rmsd_df['holo_RoG']) &
                         (filtered_rmsd_df['holo_RoG'] * MAX_ROG_PERCENTAGE > filtered_rmsd_df['apo_RoG']) &
                         (filtered_rmsd_df['holo_RoG'] * MIN_ROG_PERCENTAGE < filtered_rmsd_df['apo_RoG'])]

    return filtered_rmsd_df


def main():
    filtered_rmsd_df = load_ahojdb()
    filtered_rmsd_df = filter_utils.filter_valid_ligands(
        filtered_rmsd_df, OUTPUT_PATH)
    filter_utils.write_uniprot_ids(filtered_rmsd_df, OUTPUT_PATH)
    filter_utils.download_sequences(OUTPUT_PATH)
    filter_utils.run_shell_mmseq(OUTPUT_PATH)
    clusters = filter_utils.read_clusters(
        f'{OUTPUT_PATH}/clusterRes_cluster.tsv')
    indices = filter_utils.find_highest_cluster_pocket_rmsd_ids(
        clusters, filtered_rmsd_df)
    enhanced = filter_utils.get_other_pockets_for_each_record(
        filtered_rmsd_df, indices)
    enhanced.to_csv(DATASET_DATAFRAME)
    pockets = filter_utils.remove_duplicated_pockets(enhanced, clusters)
    filter_utils.save_dataset(pockets, OUTPUT_PATH)


if __name__ == '__main__':
    main()

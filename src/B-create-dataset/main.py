import pandas as pd
import filter_utils
import os


INPUT_PATH = '../../data/A-filter-ahojdb-v2'
INPUT_CSV = f'{INPUT_PATH}/pairs.csv'
OUTPUT_PATH = '../../data/B-create-dataset/ahoj-v2'

PREFILTERED_CSV = f'{OUTPUT_PATH}/filtered_rmsd.csv'
DATASET_DATAFRAME = f'{OUTPUT_PATH}/dataset_dataframe.csv'


def load_ahojdb():
    print('Load AHoJ-DB ...')

    # CAUTION: you might need to run a script from the /home/skrhakv/sequence-lengths-for-cryptobench to add the information about sequence length
    if not os.path.exists(PREFILTERED_CSV):
        ahojdb_df = pd.read_csv(INPUT_CSV)
        filtered_rmsd_df = ahojdb_df[(
            ahojdb_df['apo_pocket_rms'] > filter_utils.MAX_POCKET_RMSD)]
        filtered_rmsd_df.to_csv(PREFILTERED_CSV)
        del ahojdb_df
    else:
        filtered_rmsd_df = pd.read_csv(PREFILTERED_CSV)
        
    # filter AHoJ-DB according to the defined thresholds
    filtered_rmsd_df = filter_utils.get_well_defined_pairs(filtered_rmsd_df)
    filtered_rmsd_df = filtered_rmsd_df[(filtered_rmsd_df['apo_pocket_rms'] > filter_utils.MAX_POCKET_RMSD)]
    return filtered_rmsd_df


def main():
    filtered_rmsd_df = load_ahojdb()
    filtered_rmsd_df = filter_utils.filter_valid_ligands(
        filtered_rmsd_df, OUTPUT_PATH)
    filter_utils.write_uniprot_ids(filtered_rmsd_df, OUTPUT_PATH)
    filter_utils.download_sequences(OUTPUT_PATH)
    filter_utils.run_shell_mmseq(OUTPUT_PATH)
    clusters, _ = filter_utils.read_clusters(
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

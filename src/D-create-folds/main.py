import json
import sys
import csv
import shutil 

sys.path.append('../B-create-dataset')
import filter_utils

TRAIN_RATIO, TEST_RATIO = 0.8, 0.2
NUM_OF_FOLDS = 4
FOLD_RATIO = TRAIN_RATIO / NUM_OF_FOLDS

OUTPUT_PATH = '/home/vit/Projects/cryptobench/data/D-create-folds/ahoj-v2'
INPUT_PATH = '/home/vit/Projects/cryptobench/data/C-remove-holo-homomers/ahoj-v2'

def copy_shell_scripts():
    shutil.copy('../B-create-dataset/run-mmseq.sh', '.')
    shutil.copy('../B-create-dataset/download-sequences.sh', '.')

def get_uniprot_to_pdb_mapping(dataset):
    print('get UNIPROT->PDB mapping ...')

    uniprot_to_pdb = {}
    for pdb_id, holo_structures in dataset.items():
        for holo_structure in holo_structures:
            # sanity check: UNIPROT is either not present in dictionary or is already assigned to THIS PDB_ID
            assert ((str(holo_structure["uniprot_id"]) not in uniprot_to_pdb) or (
                uniprot_to_pdb[str(holo_structure["uniprot_id"])] == pdb_id))

            uniprot_to_pdb[str(holo_structure["uniprot_id"])] = pdb_id

    return uniprot_to_pdb


def write_uniprot_ids_to_file(dataset):
    print('write UNIPROT IDs to a file ... ')

    filename = f'{OUTPUT_PATH}/uniprot_ids.txt'
    with open(filename, 'w') as f:
        for pdb_id, holo_structures in dataset.items():
            for holo_structure in holo_structures:
                f.write(f'{str(holo_structure["uniprot_id"])}\n')

    # split ids which are composed of multiple chains
    with open(filename, 'r') as file:
        filedata = file.read()

    filedata = filedata.replace('-', '\n')

    with open(filename, 'w') as file:
        file.write(filedata)

    # this might create some duplicates, remove them
    with open(filename, 'r') as f:
        unique_lines = set(f.readlines())
    with open(filename, 'w') as f:
        f.writelines(unique_lines)

def get_multiple_sequences(dataset):
    # get set of multiple-sequence structures
    print('Find structures of multiple sequences ...')
    multiple_sequences = set()
    for pdb_id, holo_structures in dataset.items():
        for holo_structure in holo_structures:
            if '-' in holo_structure['uniprot_id']:
                multiple_sequences.add(holo_structure['uniprot_id'])
    return multiple_sequences

# def merge_clusters(clusters, seq_id_to_cluster, main_cluster, others, original_):
#     main_cluster_id = seq_id_to_cluster[main_cluster]
#     seq_id_to_cluster[mult_sequence] = main_cluster_id
#     clusters[main_cluster_id].append(mult_sequence)
#     for other in others:
#         other_id = seq_id_to_cluster[other]
#         extending_array = clusters[other_id]
#         clusters[main_cluster_id].extend(extending_array)
#         del clusters[other_id]
#         for i in extending_array:
#             seq_id_to_cluster[i] = main_cluster_id

def merge_clusters_with_multiple_sequence_structures(clusters, seq_id_to_cluster, multiple_sequences):
    print('Merge clusters with multiple sequence structures ...')

    for mult_sequence in multiple_sequences:
        main_cluster, others = mult_sequence.split('-')[0], mult_sequence.split('-')[1:]
        main_cluster_id = seq_id_to_cluster[main_cluster]
        seq_id_to_cluster[mult_sequence] = main_cluster_id
        clusters[main_cluster_id].append(mult_sequence)
        for other in others:
            other_id = seq_id_to_cluster[other]
            extending_array = clusters[other_id]
            clusters[main_cluster_id].extend(extending_array)
            del clusters[other_id]
            for i in extending_array:
                seq_id_to_cluster[i] = main_cluster_id
    return clusters

def assign_structures_to_folds(dataset, clusters, uniprot_to_pdb):
    print('Assign structures to folds ...')
    dataset_length = len(dataset)
    test_length = TEST_RATIO * dataset_length
    fold_length = FOLD_RATIO * dataset_length
    test = {}
    fold = -1
    train = [{} for _ in range(NUM_OF_FOLDS)]

    for key, uniprot_ids in clusters.items():
        fold = -1
        if len(test) < test_length: fold = -1
        else:
            for i in range(NUM_OF_FOLDS):
                if len(train[i]) <= fold_length:
                    fold = i
            # this might happen due to rounding error
            if fold == -1:
                fold = NUM_OF_FOLDS - 1
        for uniprot_id in uniprot_ids:
            if uniprot_id not in uniprot_to_pdb:
                continue
            pdb_id = uniprot_to_pdb[uniprot_id]
            if pdb_id == '1fzc':
                print(f'what {fold}')
            if fold == -1:
                if pdb_id not in test:
                    test[pdb_id] = dataset[pdb_id] 
            else:
                if pdb_id not in train[fold]:
                    train[fold][pdb_id] = dataset[pdb_id] 
    
    # sanity check: sum of structures in TRAIN + TEST equals the length of the whole dataset 
    fold_sum = sum([len(i) for i in train]) + len(test)
    assert fold_sum == len(dataset)

    return train, test

def merge_clusters(main_cluster_id, other_cluster_id, clusters, seq_id_to_cluster):
    extending_array = clusters[other_cluster_id]
    clusters[main_cluster_id].extend(extending_array)
    del clusters[other_cluster_id]
    for i in extending_array:
        seq_id_to_cluster[i] = main_cluster_id
    
    return clusters, seq_id_to_cluster

def check_structures_assigned_to_multiple_uniprot_ids(uniprot_to_pdb_mapping, clusters, seq_id_to_cluster):
    pdb_to_uniprot = {}
    for uniprot_id, pdb_id in uniprot_to_pdb_mapping.items():
        if pdb_id not in pdb_to_uniprot:
            pdb_to_uniprot[pdb_id] = uniprot_id
        else:
            # check that the PDB_IDs are in the same cluster:
            uniprot_id1 = uniprot_id
            uniprot_id2 = pdb_to_uniprot[pdb_id]

            cluster_id1 = seq_id_to_cluster[uniprot_id1] 
            cluster_id2 = seq_id_to_cluster[uniprot_id2]
            if cluster_id1 != cluster_id2:
                clusters, seq_id_to_cluster = merge_clusters(cluster_id1, cluster_id2, clusters, seq_id_to_cluster)

            assert seq_id_to_cluster[uniprot_id1] == seq_id_to_cluster[uniprot_id2], f'{pdb_id}, {uniprot_id1}, {uniprot_id2}' 

def save_dataset(train, test):
    print('Save dataset ...')
    with open(f'{OUTPUT_PATH}/cryptobench/folds/test.json', 'w', encoding='utf-8') as f:
        json.dump(test, f, ensure_ascii=False, indent=4)
    for i in range(NUM_OF_FOLDS):
        with open(f'{OUTPUT_PATH}/cryptobench/folds/train-fold-{i}.json', 'w', encoding='utf-8') as f:
            json.dump(train[i], f, ensure_ascii=False, indent=4)

    splits = {}
    splits['test'] = list(test.keys())
    for i in range(NUM_OF_FOLDS):
        splits[f'train-{i}'] = list(train[i].keys())
    with open(f'{OUTPUT_PATH}/cryptobench/splits.json', 'w', encoding='utf-8') as f:
        json.dump(splits, f, ensure_ascii=False, indent=4)

def main():
    with open(f'{INPUT_PATH}/dataset.json') as f:
        dataset = json.load(f)
    copy_shell_scripts()

    uniprot_to_pdb_mapping = get_uniprot_to_pdb_mapping(dataset)
    # write_uniprot_ids_to_file(dataset)
    # filter_utils.download_sequences(OUTPUT_PATH)
    # filter_utils.run_shell_mmseq(OUTPUT_PATH)
    multiple_sequences = get_multiple_sequences(dataset)
    clusters, seq_id_to_cluster = filter_utils.read_clusters(f'{OUTPUT_PATH}/clusterRes_cluster.tsv')
    clusters = merge_clusters_with_multiple_sequence_structures(clusters, seq_id_to_cluster, multiple_sequences)

    # some PDB IDs might be assigned to multiple sequences! therefore it is necessary to merge those clusters together
    check_structures_assigned_to_multiple_uniprot_ids(uniprot_to_pdb_mapping, clusters, seq_id_to_cluster)
    train, test = assign_structures_to_folds(dataset, clusters, uniprot_to_pdb_mapping)
    save_dataset(train, test)

if __name__ == '__main__':
    main()

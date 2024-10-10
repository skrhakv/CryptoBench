import subprocess
import csv
import pandas as pd
import pocket_selection_parser
import pocket
import warnings
import json
from rdkit import Chem
import os

# # #1 denotes hydrogen, #6 carbon, etc.
# ORGANIC_ATOM_QUERY = Chem.MolFromSmarts(
#     '[!$([#1,#5,#6,#7,#8,#9,#15,#16,#17,#35,#53])]')
# # check for the existence of any kind of CC bond
# ORGANIC_BOND_QUERY = Chem.MolFromSmarts('[#6]-,=,#,:[#6]')
# HAS_CH_QUERY = Chem.MolFromSmarts('[C!H0]')

cached_smiles = {}

# https://github.com/rdkit/rdkit/discussions/5314

P2RANK_ATOMS_NUM_THRESHOLD = 5

# DEFINE FILTER THRESHOLDS
MAX_POCKET_RMSD = 0.5
MIN_TM_SCORE = 0.5
MAX_POCKET_DISTANCE = 4
MIN_ROG_PERCENTAGE = 0.8
MAX_ROG_PERCENTAGE = 1.2
MIN_SEQUENCE_LENGTH = 50

def get_well_defined_pairs(df):
    return df[(df['apo_tm_score'] > MIN_TM_SCORE) &
                         (df['apo_tm_score_i'] > MIN_TM_SCORE) &
                         (df['apo_pocket_dist'] < MAX_POCKET_DISTANCE) &
                         (df['apo_RoG'] * MAX_ROG_PERCENTAGE > df['holo_RoG']) &
                         (df['apo_RoG'] * MIN_ROG_PERCENTAGE < df['holo_RoG']) &
                         (df['holo_RoG'] * MAX_ROG_PERCENTAGE > df['apo_RoG']) &
                         (df['holo_RoG'] * MIN_ROG_PERCENTAGE < df['apo_RoG']) &
                         # check that each sequence comply the MIN_SEQUENCE_LENGTH requirement
                         (df['apo_UNPovrlp_obs'] > MIN_SEQUENCE_LENGTH)]
                         # (df['apo_sequence_length'].apply(lambda sequence_lengths: all(int(i) > MIN_SEQUENCE_LENGTH for i in str(sequence_lengths).split('_'))))]

def check_ligand_atom_count(smiles):
    # TODO: switch for p2rank-based filtering:
    # hint: https://chemistry.stackexchange.com/questions/43363/how-to-calculate-hydrogen-from-smiles-string
    # TODO: test manually on Zn++, etc.
    if len(smiles.values) < 1:
        return False

    for i in str(smiles.values[0]).split(';'):
        i = i.strip()
        if i in cached_smiles.keys():
            if cached_smiles[i]:
                return True
            else:
                continue
        try:
            molecule = Chem.MolFromSmiles(i)
            atoms_count = molecule.GetNumAtoms()
        except:
            cached_smiles[i] = False
            continue

        is_valid_smiles = True
        if atoms_count < P2RANK_ATOMS_NUM_THRESHOLD:
            is_valid_smiles = False

        cached_smiles[i] = is_valid_smiles
        return is_valid_smiles

        #
        # this is the old way - look into the BioLip database to find its smiles and then use the ORGANIC_ATOM_QUERY, ORGANIC_BOND_QUERY, HAS_CH_QUERY whether
        # the ligand is 'organic' or not
        #

        # try:
        #     _is_organic = (not molecule.HasSubstructMatch(ORGANIC_ATOM_QUERY)) and molecule.HasSubstructMatch(ORGANIC_BOND_QUERY)\
        #         and molecule.HasSubstructMatch(HAS_CH_QUERY)
        #     cached_smiles[i] = _is_organic
        #     return _is_organic
        # except:
        #     cached_smiles[i] = False
        #     continue
    return False


IGNORED_GROUPS_LIST = ['HOH', 'DOD', 'WAT', 'UNK', 'ABA', 'MPD', 'GOL', 'SO4', 'PO4']

def remove_ignored_groups(df):
    """P2RANK defines a list of ignored groups: (HOH, DOD, WAT, NAG, MAN, UNK, GLC, ABA, MPD, GOL, SO4, PO4). 
    However, we would like to keep NAG, GLC, MAN. This function removes the ignored group from the dataframe except for NAG, GLC, MAN."""
    return df[~df['ligand'].isin(IGNORED_GROUPS_LIST)]


# def remove_multichain_pockets(df):
#     return df[
#         (df['apo_chains3'].str.contains('-') != True) &
#         (df['holo_chains3'].str.contains('-') != True) &
#         (df['apo_UNPs'].str.contains('-') != True) &
#         (df['holo_UNPs'].str.contains('-') != True) &
#         (df['apo_structure'].str.contains('\+') != True) &
#         (df['holo_structure'].str.contains('\+') != True)]


def filter_valid_ligands(df, path):
    print('Filter ligands ...')
    cached_smiles.clear()

    # read cached smiles
    cached_smiles_path = f'{path}/cached_smiles.csv'
    if os.path.exists(cached_smiles_path):
        with open(cached_smiles_path, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=';')
            for row in reader:
                cached_smiles[row[0].strip()] = row[1] == 'True'

    # load all ligand smiles
    valid_ligands_df = pd.read_csv(f'{path}/ligand.tsv', sep='\t')
    df[['ligand_chain', 'ligand', 'ligand_index']
       ] = df['apo_query_POI'].str.split('_', expand=True)

    # check each ligand in dataset using smiles
    df = remove_ignored_groups(df)
    filtered_df = df[
        (df['ligand'].isin(valid_ligands_df['#CCD'].values)) &
        (df['ligand'].apply(lambda ligand: check_ligand_atom_count(valid_ligands_df[valid_ligands_df['#CCD'] == ligand]['SMILES'])))]

    # save cached results for better performance next time
    with open(cached_smiles_path, 'w') as f:
        for key, value in cached_smiles.items():
            f.write(f"{key.strip()};{value}\n")

    return filtered_df


def write_uniprot_ids(filtered_rmsd_df, output_path):
    with open(f'{output_path}/uniprot_ids.txt', 'w') as f:
        for uniprot_id in filtered_rmsd_df['apo_UNPs'].unique():
            ids = str(uniprot_id).split(' ')
            for i in ids:
                f.write(f'{i}\n')


def download_sequences(path):
    """downloads sequences from UniProt"""
    print('Download sequences ...')

    subprocess.call(['sh', './download-sequences.sh', path])


def run_shell_mmseq(path, min_seq_identity=0.4):
    """Runs mmseqs2 in shell"""
    print('Run mmseqs ...')
    subprocess.call(['sh', './run-mmseq.sh', path, str(min_seq_identity)])


def read_clusters(clusters_filepath):
    print('Read clusters ...')

    seq_id_to_cluster = {}
    clusters = {}
    with open(clusters_filepath, 'r') as csvfile:
        filtered_uniprot_ids = csv.reader(csvfile, delimiter='\t')
        for row in filtered_uniprot_ids:
            if row[0] not in clusters:
                clusters[row[0]] = []
            clusters[row[0]].append(row[1])
            seq_id_to_cluster[row[1]] = row[0]

    return clusters, seq_id_to_cluster


def find_highest_cluster_pocket_rmsd_ids(clusters, df):
    print('Find highest cluster pocket rmsds ...')

    indices = []
    for cluster, cluster_sequences in clusters.items():
        indices.append(
            df[
                df['apo_UNPs'].apply(lambda s: any(x in s for x in cluster_sequences))]['apo_pocket_rms'].idxmax())
    return indices


def get_other_pockets_for_each_record(df, largest_rmsd_indices):
    print('Enhance the selected dataset with other pockets ...')

    largest_pockets_for_unique_uniprot = df.loc[largest_rmsd_indices]
    id_with_chain = list(zip(
        largest_pockets_for_unique_uniprot['apo_structure'], largest_pockets_for_unique_uniprot['apo_chains3']))

    additional_records = []
    # we want to search for pockets in single chains to enhance multi-chain structures:
    for structure, chain in id_with_chain:
        if '-' in chain:
            for single_chain in chain.split('-'):
                additional_records.append((structure, single_chain))

    id_with_chain.extend(additional_records)

    tuples_in_df = pd.MultiIndex.from_frame(
        df[['apo_structure', 'apo_chains3']])

    enhanced = df[tuples_in_df.isin(id_with_chain)]
    enhanced = enhanced[(enhanced['holo_pocket_selection'].notna())]
    return enhanced


def remove_duplicated_pockets(enhanced, clusters):
    print('Remove duplicates ...')

    grouped_df = enhanced.groupby(['apo_structure', 'apo_chains3'])

    warnings.filterwarnings("ignore")

    parser = pocket_selection_parser.pocket_selection_parser()
    pockets: dict[tuple[str, str], list[pocket.pocket]] = {}

    used_keys = set()
    
    # for each group (group=pdb_id + chain_id) parse its pocket
    for key, item in grouped_df:
        group = grouped_df.get_group(key)

        # merge multichain entries with their singlechain subsets
        if '-' in key[1]:
            for i in key[1].split('-'):
                new_key = (key[0], i)
                if new_key in grouped_df.groups.keys():
                    group = pd.concat([group, grouped_df.get_group(new_key)])
                    used_keys.add(new_key)
        
        # parse
        pdb_id, pockets_for_this_pdb_id = parser.parse_pocket_selections(group)
        pockets[key] = pockets_for_this_pdb_id

    for i in used_keys:
        del pockets[i]

    duplicates = []
    for key, value in clusters.items():
        if len(value) > 1:
            for i in range(1, len(value)):
                duplicates.append((value[i], False))

    remove_duplicates = set()
    for id, pckts in pockets.items():
        for index, item in enumerate(duplicates):
            if not item[1] and item[0] in pckts[0].uniprot_id:
                duplicates[index] = (duplicates[index][0], True)
                remove_duplicates.add(id)
    for i in remove_duplicates:
        del pockets[i]
    return pockets


def save_dataset(pockets, output_path, filename='dataset.json'):
    print('Save the dataset ...')

    FOR_PYMOL = True
    data = {}
    for apo_pdb_id, pdb_pockets in pockets.items():
        apo_pdb_id = apo_pdb_id[0]
        for p in pdb_pockets:
            if apo_pdb_id not in data:
                data[apo_pdb_id] = []
            if FOR_PYMOL:
                data[apo_pdb_id].append(
                    {'uniprot_id': p.uniprot_id.replace(' ', '-'), 'holo_pdb_id': p.holo_pdb_id, 'holo_chain': p.holo_chain, 'apo_chain': p.chain,
                     'ligand': p.ligand, 'ligand_index': p.ligand_index, 'ligand_chain': p.ligand_chain, 'pRMSD': p.pocket_rmsd,
                     'apo_pocket_selection': p.get_apo_pocket_definition(), 'holo_pocket_selection': p.get_holo_pocket_definition(),
                     'apo_pymol_selection': p.apo_selection, 'holo_pymol_selection': p.holo_selection})
            else:
                data[apo_pdb_id].append(
                    {'uniprot_id': p.uniprot_id.replace(' ', '-'), 'holo_pdb_id': p.holo_pdb_id, 'holo_chain': p.holo_chain, 'apo_chain': p.chain,
                     'ligand': p.ligand, 'ligand_index': p.ligand_index, 'ligand_chain': p.ligand_chain, 'pRMSD': p.pocket_rmsd,
                     'apo_pocket_selection': p.get_apo_pocket_definition(), 'holo_pocket_selection': p.get_holo_pocket_definition()})

    with open(f'{output_path}/{filename}', 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)

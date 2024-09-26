
import biotite.database.rcsb as rcsb
import os
from biotite.structure.io.pdbx import get_structure
import biotite.structure.io.pdbx as pdbx
import biotite
import numpy as np

CIF_FILES_PATH = '/home/vit/Projects/deeplife-project/data/cif_files'

# pdb_id = '4wmb'
# chain_id = 'D'


def get_sequence_length(pdb_id, chain_id):
    cif_file_path = rcsb.fetch(pdb_id, "cif", CIF_FILES_PATH)
    cif_file = pdbx.CIFFile.read(cif_file_path)

    structure = get_structure(cif_file, model=1)
    all_atoms = structure[(structure.chain_id == chain_id) & (biotite.structure.filter_peptide_backbone(structure))]


    chain_id_changes = (all_atoms.chain_id[1:] != all_atoms.chain_id[:-1])
    res_id_changes   = (all_atoms.res_id[1:]   != all_atoms.res_id[:-1]  )
    ins_code_changes = (all_atoms.ins_code[1:] != all_atoms.ins_code[:-1])
    res_name_changes = (all_atoms.res_name[1:] != all_atoms.res_name[:-1])
    residue_change_mask = (
        chain_id_changes |
        res_id_changes |
        ins_code_changes |
        res_name_changes
    )

    residue_starts = np.where(residue_change_mask)[0] + 1
    residue_starts = np.concatenate(([0], residue_starts))

    return len(residue_starts)




# import requests
# 
# pdb_id = '4wmb'
# chain_id = 'D'
# 
# CIF_FILES_PATH = '/home/vit/Projects/deeplife-project/data/cif_files'
# 
# def get_entity_id(pdb_id, chain_id):
#     pdb_info = requests.get(
#         f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}").json()
#     print(pdb_info)
#     for entity in pdb_info[pdb_id]:
#         if chain_id.upper() in entity['in_chains']:
#             return entity['entity_id']
#     return None
# 
# 
# def get_sequence_length(pdb_id, entity_id):
#     pdb_info = requests.get(
#         f'https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/uniprot_mapping/{pdb_id}/{entity_id}')
#     print(pdb_info[pdb_id])
# 
# entity_id = get_entity_id(pdb_id, chain_id)
# get_sequence_length(pdb_id, entity_id)
import json
from itertools import groupby
import biotite.database.rcsb as rcsb
import biotite.structure.io.pdbx as mmcif
import biotite.structure as struc

CIF_FILES_PATH = '/home/vit/Projects/deeplife-project/data/cif_files'
#INPUT_PATH = '/home/vit/Projects/cryptobench/data/I-final/cryptobench/cryptobench-dataset/dataset.json'
INPUT_PATH = '/home/vit/Projects/cryptobench/data/I-final/cryptobench/cryptobench-dataset/auxiliary-data/non-cryptic-pockets/noncryptic-pockets.json'
NON_CRYPTIC = True

def auth_to_label(pdb_id, auth_binding_residues):
    chain_id = auth_binding_residues[0].split('_')[0]

    # download and read the cif file
    cif_file_path = rcsb.fetch(pdb_id, "cif", CIF_FILES_PATH)
    cif_file = mmcif.PDBxFile.read(cif_file_path)

    # load auth
    whole_structure = mmcif.get_structure(cif_file, model=1)
    protein = whole_structure[struc.filter_amino_acids(whole_structure)]
    auth_ca_atoms = protein[(protein.atom_name == "CA") &
                            (protein.element == "C")]

    # load label
    whole_structure = mmcif.get_structure(
        cif_file, model=1, use_author_fields=False)
    protein = whole_structure[struc.filter_amino_acids(whole_structure)]
    label_ca_atoms = protein[(protein.atom_name == "CA") & (
        protein.element == "C")]

    # sanity check
    assert (len(auth_ca_atoms) == len(label_ca_atoms))

    auth_binding_residues = [i.split('_')[1] for i in auth_binding_residues]
    label_binding_residues = []

    for auth_ca, label_ca in zip(auth_ca_atoms, label_ca_atoms):
        # skip residues which are not in desired chain
        if auth_ca.chain_id != chain_id:
            continue
        
        # sanity check
        assert (auth_ca.res_name == label_ca.res_name)
        if str(auth_ca.res_id) in auth_binding_residues:
            label_binding_residues.append(
                f'{chain_id}_{label_ca.res_id}')
    
    return label_binding_residues

    


def main():
    with open(INPUT_PATH) as f:
        ds = json.load(f)

    label_ds = {}

    for apo_pdb_id, holo_structures in ds.items():
        auth_binding_residues = set()
        for holo_structure in holo_structures:
            auth_binding_residues.update(holo_structure['apo_pocket_selection'])

        # group by chain id
        auth_binding_residues = [list(g) for k, g in groupby(auth_binding_residues, key=lambda x: x.split('_')[0])]
        
        label_binding_residues = []
        for chain_group in auth_binding_residues:
            label_binding_residues.extend(auth_to_label(apo_pdb_id, chain_group))
        label_ds[apo_pdb_id] = list(label_binding_residues)

    if NON_CRYPTIC:
        with open('label_dataset.json') as f:
            ds = json.load(f)
        new_label_ds = {}
        for apo_pdb_id in label_ds.keys():
            cryptic_binding_residues = ds[apo_pdb_id]
            noncryptic_binding_residues = label_ds[apo_pdb_id]

            if len([i for i in noncryptic_binding_residues if i in cryptic_binding_residues ]) / len(noncryptic_binding_residues) > 0.75:
                # new_label_ds[apo_pdb_id] = []
                continue
            else: new_label_ds[apo_pdb_id] = [i for i in noncryptic_binding_residues if i not in cryptic_binding_residues]
        

    output_file = 'label_dataset.json' if not NON_CRYPTIC else 'label_noncryptic_dataset.json'
    with open(output_file, 'w') as f:
        json.dump(new_label_ds, f)

if __name__ == '__main__':
    main()
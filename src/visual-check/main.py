from pymol import cmd
import time
import json
import sys
import pymol

_stdouterr = sys.stdout, sys.stderr
pymol.finish_launching()
sys.stdout, sys.stderr = _stdouterr

DATASET_PATH = '/home/vit/Projects/cryptobench/data/B-create-dataset/single-chain-with-covalent-ligands/dataset.json'
CIF_FILES_PATH = '/home/vit/Projects/ahoj2-extraction/src/analysis/cif_files'

with open(DATASET_PATH) as f:
    dataset = json.load(f)

index = 0
for pdb_id, data in dataset.items():
    index += 1
    data= data[0]
    cmd.reinitialize()
    cmd.set('fetch_path', cmd.exp_path(CIF_FILES_PATH), quiet=0)
    try:
        cmd.fetch(f'{pdb_id}{data["apo_chain"]}')
        cmd.fetch(f'{data["holo_pdb_id"]}{data["holo_chain"]}')
        
    except:
        continue

    cmd.color('blue', f'{pdb_id}{data["apo_chain"]}')
    cmd.show('sticks', data['apo_pymol_selection'])

    cmd.color('forest', f'{data["holo_pdb_id"]}{data["holo_chain"]}')
    cmd.show('sticks', data['holo_pymol_selection'])
    try:
        cmd.color('aquamarine', data['apo_pymol_selection'])
        cmd.color('yellow', data['holo_pymol_selection'])
    except:
        continue

    results1 = cmd.align(f'{pdb_id}{data["apo_chain"]}', f'{data["holo_pdb_id"]}{data["holo_chain"]}')

    print(pdb_id + data["apo_chain"])
    cmd.zoom(f'{data['apo_pymol_selection']}')

    time.sleep(10)




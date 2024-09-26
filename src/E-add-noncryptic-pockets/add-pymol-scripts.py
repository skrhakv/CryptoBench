import json
from itertools import groupby

OUTPUT_PATH = '/home/vit/Projects/cryptobench/data/E-add-noncryptic-pockets/ahoj-v2/cryptobench'
CRYPTIC_DATASET_PATH = f'{OUTPUT_PATH}/dataset.json'
NONCRYPTIC_DATASET_PATH = f'{OUTPUT_PATH}/auxiliary-data/non-cryptic-pockets/noncryptic-pockets.json'
PYMOL_SCRIPTS_PATH = f'{OUTPUT_PATH}/auxiliary-data/pymol-scripts'

def create_id(apo_pdb_id, record):
    return f'{apo_pdb_id}{record["apo_chain"]}_{record["holo_pdb_id"]}{record["holo_chain"]}'

class Pocket:
    def __init__(self, apo_pdb_id, cryptic_record) -> None:
        self.id = create_id(apo_pdb_id, cryptic_record)
        [self.apo_id, self.holo_id] = self.id.split('_')
        self.apo_pdb_id = self.apo_id[:4]
        self.holo_pdb_id = self.holo_id[:4]
        self.apo_chains = self.apo_id[4:].split('-')
        self.holo_chains = self.holo_id[4:].split('-')
        self.cryptic_apo_selection = set(cryptic_record['apo_pocket_selection'])
        self.cryptic_holo_selection = set(cryptic_record['holo_pocket_selection'])
    
    def add_noncryptic_pocket(self, record):
        if 'noncryptic_apo_selection' in dir(self) and 'noncryptic_holo_selection' in dir(self):
            self.extend_noncryptic_pocket(record)
        else:
            self.noncryptic_apo_selection = set(record['apo_pocket_selection'])
            self.noncryptic_holo_selection = set(record['holo_pocket_selection'])

    def extend_cryptic_pocket(self, record):
        assert self.cryptic_apo_selection is not None
        assert self.cryptic_holo_selection is not None
        self.cryptic_apo_selection.update(record['apo_pocket_selection'])
        self.cryptic_holo_selection.update(record['holo_pocket_selection'])

    def extend_noncryptic_pocket(self, record):
        assert self.noncryptic_apo_selection is not None
        assert self.noncryptic_holo_selection is not None
        self.noncryptic_apo_selection.update(record['apo_pocket_selection'])
        self.noncryptic_holo_selection.update(record['holo_pocket_selection'])

    def get_pymol_selection(self, pdb_id, selection):
        groups = [list(group) for _, group in groupby(selection, key=lambda x: x.split('_')[0])]
        chain_selections = [] 
        for group in groups:
            chain_id = group[0].split('_')[0]
            residues = [i.split('_')[1] for i in group]
            selection = f'( chain {chain_id} and resi {"+".join(residues)} )'
            chain_selections.append(selection)
        return f'{pdb_id} and ( {" or ".join(chain_selections)} )'


def main():
    with open(CRYPTIC_DATASET_PATH) as f:
        cryptic_dataset = json.load(f)

    with open(NONCRYPTIC_DATASET_PATH) as f:
        noncryptic_dataset = json.load(f)

    pockets: {str, Pocket} = {}
    for apo_pdb_id, holo_structures in cryptic_dataset.items():
        for holo_structure in holo_structures:
            id = create_id(apo_pdb_id, holo_structure)
            if id in pockets:
                pockets[id].extend_cryptic_pocket(holo_structure)
            else: 
                pockets[id] = Pocket(apo_pdb_id, holo_structure)
    
    for apo_pdb_id, holo_structures in noncryptic_dataset.items():
        for holo_structure in holo_structures:
            id = create_id(apo_pdb_id, holo_structure)
            if id in pockets:
                pockets[id].add_noncryptic_pocket(holo_structure)
                
    for id, pocket in pockets.items():
        with open(f'{PYMOL_SCRIPTS_PATH}/{pocket.id}.pml', 'w') as pymol_script:

            pymol_script.write(f'reinitialize\n')

            if '-' in pocket.apo_id or '-' in pocket.holo_id:
                pymol_script.write(f'fetch {pocket.apo_pdb_id}\n')
                pymol_script.write(f'fetch {pocket.holo_pdb_id}\n')

            else:
                pymol_script.write(f'fetch {pocket.apo_id}\n')
                pymol_script.write(f'fetch {pocket.holo_id}\n')

            pymol_script.write(f'color blue, {pocket.apo_pdb_id}\n')
            if 'noncryptic_apo_selection' in dir(pocket):
                pymol_script.write(f'color white, {pocket.get_pymol_selection(pocket.apo_pdb_id, pocket.noncryptic_apo_selection)}\n')
            pymol_script.write(f'color aquamarine, {pocket.get_pymol_selection(pocket.apo_pdb_id, pocket.cryptic_apo_selection)}\n')
            pymol_script.write(f'show sticks, {pocket.get_pymol_selection(pocket.apo_pdb_id, pocket.cryptic_apo_selection)}\n')

            pymol_script.write(f'color green, {pocket.holo_pdb_id}\n')
            if 'noncryptic_holo_selection' in dir(pocket):
                pymol_script.write(f'color white, {pocket.get_pymol_selection(pocket.holo_pdb_id, pocket.noncryptic_holo_selection)}\n')
            pymol_script.write(f'color yellow, {pocket.get_pymol_selection(pocket.holo_pdb_id, pocket.cryptic_holo_selection)}\n')
            pymol_script.write(f'show sticks, {pocket.get_pymol_selection(pocket.holo_pdb_id, pocket.cryptic_holo_selection)}\n')
            
            if '-' in pocket.apo_id or '-' in pocket.holo_id:
                apo_chain_selection = " or ".join([f"chain {i}" for i in pocket.apo_chains])
                holo_chain_selection = " or ".join([f"chain {i}" for i in pocket.holo_chains])
                pymol_script.write(f'align {pocket.apo_pdb_id} and ( {apo_chain_selection} ),{pocket.holo_pdb_id} and ( {holo_chain_selection} ) \n')

            else:
                pymol_script.write(f'align {pocket.apo_pdb_id} ,{pocket.holo_pdb_id}\n')
            pymol_script.write(f'zoom {pocket.apo_pdb_id}\n')
            pymol_script.write("print 'COLORING:\\nAPO structure - blue\\nAPO cryptic region - aquamarine\\nHOLO structure - green\\nHOLO cryptic region - yellow\\nNONCRYPTIC pockets - white'")
if __name__ == '__main__':
    main()

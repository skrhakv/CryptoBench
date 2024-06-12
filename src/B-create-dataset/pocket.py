# TODO: refactor

class pocket:
    def __init__(self, pdb_id, ligand, pocket_rmsd, holo_pdb_id, chain, holo_chain, uniprot_id, apo_selection, holo_selection, ligand_index, ligand_chain):
        self.pdb_id: str = pdb_id
        self.chain: str = chain
        self.ligand: str = ligand
        self.apo_binding_residues: {str, [str]} = {}
        self.holo_binding_residues: {str, [str]} = {}
        self.pocket_rmsd = pocket_rmsd
        self.holo_pdb_id: str = holo_pdb_id
        self.holo_chain: str = holo_chain
        self.uniprot_id: str = uniprot_id
        self.apo_selection: str = apo_selection
        self.holo_selection: str = holo_selection
        self.ligand_index: str = ligand_index
        self.ligand_chain: str = ligand_chain
    def add_chain_holo_binding_residues(self, chain: str, binding_residues: [str]):
        if not chain in self.holo_binding_residues:
            self.holo_binding_residues[chain] = binding_residues
        else: self.holo_binding_residues[chain].extend(binding_residues)

    def add_chain_apo_binding_residues(self, chain: str, binding_residues: [str]):
        if not chain in self.apo_binding_residues:
            self.apo_binding_residues[chain] = binding_residues
        else: self.apo_binding_residues[chain].extend(binding_residues)

    def get_apo_pocket_definition(self) -> [str]:
        pckt_def: [str] = []
        for chain, definition in self.apo_binding_residues.items():
            pckt_def.extend([f'{chain}_{x}' for x in definition])
        return pckt_def
    
    def get_holo_pocket_definition(self) -> [str]:
        pckt_def: [str] = []
        for chain, definition in self.holo_binding_residues.items():
            pckt_def.extend([f'{chain}_{x}' for x in definition])
        return pckt_def

    def __str__(self):
        return f'{self.pdb_id} {self.ligand} {self.binding_residues}'

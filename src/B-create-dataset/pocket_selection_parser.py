import pandas as pd
from pocket import pocket

# TODO: refactor

class pocket_selection_parser:
    def parse_pocket_selections(self, apo_holo_pairs: pd.DataFrame) -> (str, [pocket]):
        self.__get_ligands(apo_holo_pairs)
        # selections = apo_holo_pairs['apo_pocket_selection']
        selections = list(zip(
            apo_holo_pairs['ligand'],
            apo_holo_pairs['apo_pocket_selection'],
            apo_holo_pairs['apo_pocket_rms'],
            apo_holo_pairs['holo_structure'],
            apo_holo_pairs['apo_chains3'],
            apo_holo_pairs['holo_chains3'],
            apo_holo_pairs['apo_UNPs'],
            apo_holo_pairs['holo_pocket_selection'],
            apo_holo_pairs['ligand_index'],
            apo_holo_pairs['ligand_chain']))
        pockets: [pocket] = []
        for ligand, selection, pocket_rmsd, holo_pdb_id, chain, holo_chain, uniprot_id, holo_selection, ligand_index, ligand_chain in selections:
            pockets.append(self.__get_pocket_definition(ligand, selection, int(
                pocket_rmsd), holo_pdb_id, chain, holo_chain, uniprot_id, holo_selection, ligand_index, ligand_chain))
        return pockets[0].pdb_id, pockets

    def __get_pocket_definition(self, ligand, apo_selection: str, pocket_rmsd, holo_pdb_id, chain, holo_chain, uniprot_id, holo_selection, ligand_index, ligand_chain) -> pocket:
        apo_split_by_and = apo_selection.split(' and ')
        apo_pdb_id = apo_split_by_and[0].strip()

        holo_split_by_and = holo_selection.split(' and ')
        parsed_holo_pdb_id = holo_split_by_and[0].strip()
        
        assert parsed_holo_pdb_id == holo_pdb_id

        p = pocket(apo_pdb_id, ligand, pocket_rmsd, holo_pdb_id,
                   chain, holo_chain, uniprot_id, apo_selection, holo_selection, ligand_index, ligand_chain)
        
        apo_pocket_selection = ' and '.join(apo_split_by_and[1:])
        apo_single_chain_selections = apo_pocket_selection.split(' or ')
        for single_chain_selection in apo_single_chain_selections:
            self.__get_chain_definition(p, single_chain_selection, 'apo')

        holo_pocket_selection = ' and '.join(holo_split_by_and[1:])
        holo_single_chain_selections = holo_pocket_selection.split(' or ')
        for single_chain_selection in holo_single_chain_selections:
            self.__get_chain_definition(p, single_chain_selection, 'holo')
        return p

    def __get_chain_definition(self, p: pocket, selection: str, apo_holo):
        selection = selection.replace('(', '').replace(')', '').strip()
        chain_definition = selection.split(' and ')[0].strip()
        pocket_definition = selection.split(' and ')[1].strip()
        chain = chain_definition.split(' ')[1]
        pocket_selection = [s for s in pocket_definition.split(' ')[
            1].split('+')]
        if apo_holo == 'apo':
            p.add_chain_apo_binding_residues(chain, pocket_selection)
        else: p.add_chain_holo_binding_residues(chain, pocket_selection)

    def __get_ligands(self, apo_holo_pairs):
        apo_holo_pairs[['ligand_chain', 'ligand', 'ligand_index']
                       ] = apo_holo_pairs['apo_query_POI'].str.split('_', expand=True)

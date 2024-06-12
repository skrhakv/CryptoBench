from pymol import cmd

structure_set = set()


def get_pocket_residues(pocket_selection) -> [str]:
    pocket_residues = []
    previous_atom = ('', '')
    for atom in [(a.resn, a.resi) for a in cmd.get_model(selection=pocket_selection).atom]:
        if previous_atom != atom:
            pocket_residues.append(atom[0])
    return pocket_residues


def check_pocket_similarity(apo_pocket, holo_pocket) -> bool:
    apo_pocket_residues = get_pocket_residues(apo_pocket)
    holo_pocket_residues = get_pocket_residues(holo_pocket)

    if len(apo_pocket_residues) != len(holo_pocket_residues):
        return False

    for apo_pocket_residue, holo_pocket_residue in zip(apo_pocket_residues, holo_pocket_residues):
        if apo_pocket_residue != holo_pocket_residue:
            return False
    return True


def get_pocket_rmsd(apo_pocket, holo_pocket) -> int:
    return cmd.align(apo_pocket, holo_pocket, cycles=0)[0]


def load_pair(structure1, structure2) -> bool:
    if len(structure_set) > 50:
        reinitialize_pymol()
    if structure1 == structure2:
        return False
    if structure1 not in structure_set:
        cmd.fetch(structure1, path=f'tmp-data/')
        structure_set.add(structure1)
    if structure2 not in structure_set:
        cmd.fetch(structure2, path=f'tmp-data/')
        structure_set.add(structure2)

    return True


def reinitialize_pymol():
    structure_set.clear()
    cmd.reinitialize()


def get_pocket_selection(pocket_selections, structure_info) -> str | None:
    structure_id = structure_info.structure
    pocket_id = structure_info.pocket[1:]
    if pocket_selections[pocket_selections[0] == f'{structure_id}.pocket_{pocket_id}'][1].empty:
        if not pocket_selections[pocket_selections[0] == f'{structure_id}.{structure_info.query_POI}'][1].empty:
            return pocket_selections[pocket_selections[0] == f'{structure_id}.{structure_info.query_POI}'][1].values[0]
        else:
            return None
    else:
        return pocket_selections[pocket_selections[0] == f'{structure_id}.pocket_{pocket_id}'][1].values[0]

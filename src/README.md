# CryptoBench
This `README.md` file describes the dataset for **CryptoBench**, a dataset and benchmark of cryptic binding sites.

## Files
- `dataset.json`: The whole CryptoBench dataset
- `folds` folder: the dataset is split into train and test subset. The train subset is further split into four folds.
- `auxiliary-data/non_cryptic_pockets/noncryptic_pockets.json`: contains additional non-cryptic pockets for each `apo` structure in the dataset. 
- `auxiliary-data/non_cryptic_pockets/pymol-scripts` folder: contains PyMOL visualization scripts for each `apo-holo` pair in the dataset. 
- `splits.json`: the list of `apo` structures in each subset.

## Format
The following files all adhere to the same format: 
- `dataset.json`
- the files in `folds` folder
- `auxiliary-data/non_cryptic_pockets/noncryptic_pockets.json`

The dataset consists of a dictionary where each key is an apo-pdb-id, and the corresponding value is an array of holo structures. This design allows a single apo structure to be associated with multiple holo structures.

The first entry looks like this:
```json
    "1a4u": [
        {
            "uniprot_id": "P10807",
            "holo_pdb_id": "1b14",
            "holo_chain": "A",
            "apo_chain": "B",
            "ligand": "NAD",
            "ligand_index": "255",
            "ligand_chain": "A",
            "apo_pocket_selection": [
                "B_12",
                "B_14",
                // other residues
                "B_189"
            ],
            "holo_pocket_selection": [
                "A_12",
                "A_14",
                // other residues
                "A_189"
            ],
            "apo_pymol_selection": "1a4u and ( (chain B and resi 12+14+15+16+17+18+37+38+62+63+64+65+91+92+93+102+106+136+137+138+151+155+181+182+183+184+186+187+188+189) )",
            "holo_pymol_selection": "1b14 and ( (chain A and resi 12+14+15+16+17+18+37+38+62+63+64+65+91+92+93+102+106+136+137+138+151+155+181+182+183+184+186+187+188+189) )",
            "pRMSD": 2.29,
            "is_main_holo_structure": false

        },
        // other holo structures associated with 1a4u

    ],
    // other apo structures
```
This entry describes an `apo` structure `1a4uB` paired with a `holo` structure `1b14A`. The `holo` structure contains a ligand `NAD` in chain `A`, index `255`, and is associated with a UNIPROT ID `P10807`. The pocket selections for the `apo` and `holo` structures are defined by `apo_pocket_selection` and `holo_pocket_selection`, respectively. For example, `B_12` indicates auth_asym_id `B` and auth_seq_id `12`. 

The pocket RMSD (**pRMSD**) equals to 2.29, however, it is not the highest pRMSD in the context of `1a4uB` (`is_main_holo_structure` is false, therefore, there is another holo structure for `1a4uB` with higher pRMSD).

## Processing the dataset
We encourage users to use the `mmCIF` format when downloading the structures, as some `PDB` files might not be available.
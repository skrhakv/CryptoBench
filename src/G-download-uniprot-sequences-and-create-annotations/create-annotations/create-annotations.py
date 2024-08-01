import json

INPUT_PATH = '/home/vit/Projects/cryptobench/data/E-add-noncryptic-pockets/ahoj-v2/cryptobench/folds'
OUTPUT_PATH = '/home/vit/Projects/cryptobench/data/G-download-uniprot-sequences-and-create-annotations/ahoj-v2/create-annotations'

DATASET_FOLDS = ['test.json',
                'train-fold-0.json',
                'train-fold-1.json',
                'train-fold-2.json',
                'train-fold-3.json']
def main():
    for fold_filename in DATASET_FOLDS:
        with open(f'{INPUT_PATH}/{fold_filename}') as f:
            dataset = json.load(f)

        results = []
        for apo_pdb_id, holo_structures in dataset.items():
            the_most_chains = holo_structures[0]['apo_chain']
            uniprot_id = holo_structures[0]['uniprot_id']
            binding_residues = set()
            for holo_structure in holo_structures:
                # check if it is multichain, possibly update the chain
                if holo_structure['apo_chain'].count('-') > the_most_chains.count('-'):
                    the_most_chains = holo_structure['apo_chain']
                binding_residues.update(holo_structure['apo_pocket_selection'])
            results.append(f'{apo_pdb_id};{the_most_chains};{uniprot_id};{" ".join(binding_residues)};UNKNOWN')

        with open(f'{OUTPUT_PATH}/{fold_filename.split(".")[0]}.csv', 'w') as f:
            for row in results:
                f.write(f'{row}\n')


if __name__ == '__main__':
    main()
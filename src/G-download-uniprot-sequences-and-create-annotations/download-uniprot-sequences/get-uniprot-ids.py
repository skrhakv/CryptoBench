import json
import sys

dataset_path = sys.argv[1]
with open(dataset_path) as f:
    dataset = json.load(f)
with open('/home/vit/Projects/cryptobench/data/G-download-uniprot-sequences-and-create-annotations/ahoj-v2/uniprot_ids.txt', 'w') as f:
    for pdb_id, data in dataset.items():
        for i in str(data[0]["uniprot_id"]).split('-'):
            f.write(f'{i}\n')
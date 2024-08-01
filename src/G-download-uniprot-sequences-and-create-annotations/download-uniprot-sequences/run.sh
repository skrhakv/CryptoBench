#!/usr/bin/env bash

OUTPUT_PATH=/home/vit/Projects/cryptobench/data/G-download-uniprot-sequences-and-create-annotations/ahoj-v2
INPUT_PATH=/home/vit/Projects/cryptobench/data/E-add-noncryptic-pockets/ahoj-v2/cryptobench

# cleanup
rm $OUTPUT_PATH/uniprot_download_error.txt
rm -rf $OUTPUT_PATH/sequences
mkdir $OUTPUT_PATH/sequences
rm $OUTPUT_PATH/sequences.fasta

# execute
python3 get-uniprot-ids.py $INPUT_PATH/dataset.json
for acc in `cat $OUTPUT_PATH/uniprot_ids.txt` ; do curl -s "https://rest.uniprot.org/uniprotkb/$acc.fasta" ; done > $OUTPUT_PATH/sequences.fasta 2>$OUTPUT_PATH/uniprot_download_error.txt
# curl  -X  GET  -i  -F "fasta=@sequences.fasta" http://proteinspector.projekty.ms.mff.cuni.cz:42013?embedder=bert --output tmp.zip 2>embedding_download_error.txt
python3 fasta-to-txt.py
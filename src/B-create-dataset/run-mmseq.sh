#!/bin/bash

OUTPUT_PATH="$1"

rm -rf $OUTPUT_PATH/tmp
rm $OUTPUT_PATH/clusterRes_*
rm $OUTPUT_PATH/out.txt
rm $OUTPUT_PATH/err.txt
/home/vit/Public/mmseqs/bin/mmseqs easy-cluster $OUTPUT_PATH/uniprot_seqs.fasta $OUTPUT_PATH/clusterRes $OUTPUT_PATH/tmp --min-seq-id 0.4 > $OUTPUT_PATH/out.txt 2> $OUTPUT_PATH/err.txt

# How to run MMseqs
This bash snippet can be used to run the mmseqs to cluster based on at least 40% similarity:
```
/home/vit/Public/mmseqs/bin/mmseqs easy-cluster concatenated-sequences.fasta clusterRes tmp --min-seq-id 0.4 > out.txt 2> err.txt
```
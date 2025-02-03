# CryptoBench: A comprehensive dataset of protein-ligand cryptic binding sites
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/skrhakv/CryptoBench/blob/master/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtae745-blue.svg)](https://doi.org/10.1093/bioinformatics/btae745)
[![PubMed](https://img.shields.io/badge/PubMed-39693053-blue.svg)](https://pubmed.ncbi.nlm.nih.gov/39693053/)
[![Dataset](https://img.shields.io/badge/Dataset%20on%20OSF-10.17605%2FOSF.IO%2FPZ4A9%20-blue.svg)](https://osf.io/pz4a9/)


Protein cryptic binding sites are sites that are spatially malformed or inaccessible in their unbound state but become visible through some external factor, such as ligand binding (see the example below). Identifying these sites is important in many applications from bioengineering to drug discovery. CryptoBench is a large-scale dataset designed to aid in the development and evaluation of new cryptic binding site prediction methods.

Illustration of the unbound state of Cobyrinic acid a,c diamide synthase (PDB ID: [4PFS](https://www.rcsb.org/structure/4pfs)), with obscured binding site. The ligand has been artificially added to highlight that it does not fit into the pocket in this state.  |  Illustration of the bound state of the same protein (PDB ID: [5IF9](https://www.rcsb.org/structure/5IF9)) shows that there actually exists a (cryptic) binding site.
:-------------------------:|:-------------------------:
![](https://github.com/skrhakv/CryptoBench/blob/master/img/4pfs.png?raw=true)  |  ![](https://github.com/skrhakv/CryptoBench/blob/master/img/5if9.png?raw=true)


## About
CryptoBench contains [over 1,000 structures](https://academic.oup.com/view-large/500101074) making it  substantially larger than any dataset available before. It can be used for training novel cryptic binding site prediction methods as it was demonstrated by training protein language model-based baseline method [within the CryptoBench manuscript](https://academic.oup.com/view-large/figure/500101098/btae745f4.tif).

The complete CryptoBench dataset, including train-test splits, CIF files, and PyMOL visualization scripts, is available on the [OSF framework](https://osf.io/pz4a9/). 

## Tutorial
To facilitate working with CryptoBench, we offer a [`tutorial/tutorial.ipynb` notebook](https://github.com/skrhakv/CryptoBench/blob/master/tutorial/tutorial.ipynb). This tutorial provides step-by-step guidance for parsing, handling train-test splits, and visualizing data within the dataset.

## Overview
1. For details on the dataset construction process and potential reproduction purposes, please refer to the `src/README.md`.
2. A framework from [this repository](https://github.com/skrhakv/apolo/tree/cryptobench-v2) was used to train the benchmark method.
3. Since the original PocketMiner code required minor adjustments to work, the forked PocketMiner repository, along with steps on how PocketMiner was evaluated on the CryptoBench test set, can be found [here](https://github.com/skrhakv/gvp).

## Benchmark method results
In the [CryptoBench study](https://academic.oup.com/bioinformatics/article/41/1/btae745/7927823), we evaluated the performance of three methods on the CryptoBench test set: the newly developed **benchmark method (pLM-NN)**, **PocketMiner**, and **P2Rank**.

| Method      | Dataset             | AUC  | AUPRC | ACC  | FPR  | TPR  | MCC  | F1 Score |
| ----------- | ------------------- | ---- | ----- | ---- | ---- | ---- | ---- | -------- |
| pLM-NN      | CB-full [^1]        | 0.86 | 0.36  | 0.93 | 0.05 | 0.48 | 0.39 | 0.92     |
| pLM-NN      | CB-PM [^2]          | 0.88 | 0.43  | 0.93 | 0.04 | 0.52 | 0.44 | 0.93     |
| PocketMiner | CB-PM               | 0.76 | 0.19  | 0.82 | 0.16 | 0.51 | 0.22 | 0.78     |
| pLM-NN      | CB-P2RANK-apo [^3]  | 0.88 | 0.42  | 0.93 | 0.04 | 0.51 | 0.43 | 0.93     |
| P2RANK      | CB-P2RANK-apo       | 0.81 | 0.21  | 0.85 | 0.14 | 0.62 | 0.27 | 0.81     |
| P2RANK      | CB-P2RANK-holo [^4] | 0.89 | 0.34  | 0.85 | 0.15 | 0.84 | 0.38 | 0.81     |

[^1]: CB-full denotes the whole CryptoBench test set.
[^2]: CB-PM denotes the subset on which PocketMiner was evaluated 
[^3]: CB-P2RANK-apo denotes the subset on which P2Rank was evaluated
[^4]: CB-P2RANK-holo denotes the holo counterparts of the CB-P2RANK-apo subset.

If you would like to evaluate your method using this dataset or compare your predictions with the benchmark, feel free to reach out! You can contact us via [GitHub Issues](https://github.com/skrhakv/CryptoBench/issues) or by email, which can be found in the paper.

## How to cite:
If you use CryptoBench, please cite [the paper](https://academic.oup.com/bioinformatics/article/41/1/btae745/7927823):

- *Vít Škrhák, Marian Novotný, Christos P Feidakis, Radoslav Krivák, David Hoksza, CryptoBench: cryptic protein–ligand binding sites dataset and benchmark, Bioinformatics, Volume 41, Issue 1, January 2025, btae745, [https://doi.org/10.1093/bioinformatics/btae745](https://doi.org/10.1093/bioinformatics/btae745)*


or, if you prefer the `BibTeX` format:

```
@article{skrhak2024cryptobench,
    author = {Škrhák, Vít and Novotný, Marian and Feidakis, Christos P and Krivák, Radoslav and Hoksza, David},
    title = {CryptoBench: Cryptic protein-ligand binding sites dataset and benchmark},
    journal = {Bioinformatics},
    pages = {btae745},
    year = {2024},
    month = {12},
    issn = {1367-4811},
    doi = {10.1093/bioinformatics/btae745},
    url = {https://doi.org/10.1093/bioinformatics/btae745},
    eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btae745/61228599/btae745.pdf},
}
```

## Contact us
If you have any questions regarding the usage of the dataset or its assembly, comparing your method against the benchmark, or if you have any suggestions, please feel free to contact us by raising [an issue!](https://github.com/skrhakv/CryptoBench/issues)

## License
This source code is licensed under the [MIT license](https://github.com/skrhakv/CryptoBench/blob/master/LICENSE).


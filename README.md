# CryptoBench: A comprehensive dataset of protein-ligand cryptic binding sites
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/skrhakv/CryptoBench/blob/master/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtae745-blue.svg)](https://doi.org/10.1093/bioinformatics/btae745)
[![PubMed](https://img.shields.io/badge/PubMed-39693053-blue.svg)](https://pubmed.ncbi.nlm.nih.gov/39693053/)
[![Dataset](https://img.shields.io/badge/Dataset%20on%20OSF-10.17605%2FOSF.IO%2FPZ4A9%20-blue.svg)](https://osf.io/pz4a9/)


Cryptic binding sites are protein pockets that are malformed or inaccessible in their unbound state but become visible upon ligand binding. Identifying these sites is crucial for drug discovery and biomedical engineering. CryptoBench provides a large-scale dataset designed to aid in the development and evaluation of cryptic binding site prediction methods.

Illustration of the unbound state of Cobyrinic acid a,c diamide synthase (PDB ID: 4PFS), showing that the binding site is not apparent. The ligand has been artificially added to highlight that it does not fit into the pocket in this state.  |  Illustration of the bound state (PDB ID: 5IF9) reveals the cryptic binding site after ligand binding.
:-------------------------:|:-------------------------:
![](https://github.com/skrhakv/CryptoBench/blob/master/img/4pfs.png?raw=true)  |  ![](https://github.com/skrhakv/CryptoBench/blob/master/img/5if9.png?raw=true)


## About
Welcome to CryptoBench, a novel and extensive dataset tailored for cryptic binding site prediction tasks. With [over 1,000 structures](https://academic.oup.com/view-large/500101074), CryptoBench is substantially larger than any previously available dataset. It can be used for training superior methods for cryptic binding site prediction, as it was already demonstrated by training a benchmark method [within the study](https://academic.oup.com/view-large/figure/500101098/btae745f4.tif).

The complete dataset, including train-test splits, CIF files, and PyMOL visualization scripts, is available on the [OSF framework](https://osf.io/pz4a9/). 

## Tutorial
To facilitate working with CryptoBench, we offer a [`tutorial/tutorial.ipynb` notebook](https://github.com/skrhakv/CryptoBench/blob/master/tutorial/tutorial.ipynb). This tutorial provides step-by-step guidance for parsing, handling train-test splits, and visualizing data within the dataset.

## Overview
1. For details on the dataset construction process and potential reproduction purposes, please refer to the `src/README.md`.
2. A framework from [this repository](https://github.com/skrhakv/apolo/tree/cryptobench-v2) was used to train the benchmark method.

## How to cite:
If you find our dataset useful, please cite [our paper](https://academic.oup.com/bioinformatics/article/41/1/btae745/7927823):


- *Vít Škrhák, Marian Novotný, Christos P Feidakis, Radoslav Krivák, David Hoksza, CryptoBench: cryptic protein–ligand binding sites dataset and benchmark, Bioinformatics, Volume 41, Issue 1, January 2025, btae745, [https://doi.org/10.1093/bioinformatics/btae745](https://doi.org/10.1093/bioinformatics/btae745)*


or, if you prefer the `BibTex` format:

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

## License
This source code is licensed under the [MIT license](https://github.com/skrhakv/CryptoBench/blob/master/LICENSE).


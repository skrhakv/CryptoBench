# Source code overview
This README file provides an overview of the pipeline used to extract the CryptoBench dataset. The extraction process is divided into several steps, allowing for easy rollback if necessary. The steps are listed alphabetically (e.g., `A-filter-ahojdb`, `B-create-dataset`, etc.).

## Download AHoJ-DB
The CryptoBench dataset is built upon AHoJ-DB, a database containing apo-holo relationships derived from individual protein-ligand interactions. You can download the database using [this link](https://owncloud.cesnet.cz/index.php/s/ED53KZc6CwsARm4). The downloaded file includes three directories, each representing a different subset of the database.

## Extraction steps
This section outlines each extraction step individually.

### A-filter-ahojdb
This step processes the batches from AHoJ-DB, extracting all apo-holo pairs, and filters them based on the following criteria:
- Resolution better than 2.5 Å,
- Matching binding residues in both the APO and HOLO states.

The filtered data is then saved as a CSV file, which can be further processed using `pandas` in the following steps.

To run this step, use the `A-filter-ahojdb/run.sh` bash script. Since there are three subsets (directories) in the downloaded AHoJ-DB file, you’ll need to run the script separately for each subset. Make sure to adjust the following parameter in the `run.sh` script to point to the particular subset:
```shell
readonly INPUT_PATH=/PATH-TO-MY-AHOJ-DB-DIRECTORY/ahoj-db/storage3/storage/praha1/home/davidhoksza/projects/ahoj/output4
```
### B-create-dataset
The second step creates the CryptoBench dataset, and it can be run by executing the `B-create-dataset/main.py` Python script. 

#### Geometric quality assurance filtering
First, to ensure the quality of the apo-holo pairs, a geometric quality assurance filter is applied with the following criteria:
1. **TM-score**: Only pairs with a `TM-score < 0.5` are excluded to ensure structural similarity,
2. **Pocket Center Distance**: The distance between the center of the binding pockets in the apo and holo states must not exceed **4 Å**. This is to filter out domain swaps and large intrachain motions, 
3. **Compactness**: A **maximum 20% change** in the radius of gyration from the holo state is allowed,
4. **Size of the mutual apo-holo region**: At least **50 observed residues** from the protein must overlap with its UniProt sequence.

#### Crypticity criterion
Once the data has been cleaned, the **crypticity criterion** is applied using **pocket RMSD (pRMSD)**. Pairs with `pRMSD > 2 Å` are considered cryptic and retained for further analysis.

#### Ligand filtering
Ligands with fewer than 5 atoms are excluded. Additionally, any ligands belonging to the following PDB groups are ignored: HOH, DOD, WAT, UNK, ABA, MPD, GOL, SO4, PO4.

#### Sequence similarity clustering and selection of a representative
Unique UniProt sequences are clustered based on ** maximum 40% sequence similarity**. From each cluster, one representative apo-holo pair is selected based on the maximal pocket RMSD, which ensures that pairs with the most significant structural changes are included.

#### Searching for additional cryptic pockets
Since one apo structure may contain multiple cryptic pockets or a single cryptic pocket may bind another type of ligand, a second search is performed within the output of the pocket RMSD filtering phase to identify these additional pockets. These newly identified pockets are included in the dataset, so each apo structure may be paired with multiple holo structures.

### C-remove-holo-homomers
In the previous step, the search for additional cryptic pockets may have introduced some duplicate homomeric holo structures in the dataset. While these duplicates do not affect the dataset’s reliability, they do add unnecessary redundancy. To address this, the script `C-remove-holo-homomers/main.py` identifies and removes any redundant holo homomers.

### D-create-folds
To prevent information leakage between the test and training subsets in the CryptoBench dataset, we assigned APO structures to each subset based on sequence similarity. First, we clustered the CryptoBench APO structures with a **maximum of 10% sequence** similarity. These clusters were then strategically combined to form a test subset and a 4-fold training subset. This approach ensures that no subset shares more than 10% sequence similarity with any other, preserving data integrity across folds.

To execute this step, run **D-create-folds/main.py**.

### E-add-noncryptic-pockets
In addition to cryptic binding sites, some APO structures may also contain non-cryptic binding sites. The script `E-add-noncryptic-pockets/main.py` generates a separate `noncryptic-pockets.json` file to store these non-cryptic binding sites, which were filtered out during earlier stages as non-cryptic.

Additionally, `add-pymol-scripts.py` is provided to create PyMOL scripts for visualizing each apo-holo pair in the dataset.

### F-statistics
This folder contains Jupyter notebooks for calculating statistics, tables, and generating graphs used in the manuscript:
1. `compare-metrics-figure.ipynb` – Compares various pocket-oriented metrics between the PocketMiner test set and the AHoJ-DB database to highlight differences in metric distributions between a cryptic dataset and a general dataset (mix of cryptic and rigid binding sites),
2. `pocketminer-cryptosite-stats.ipynb` – Provides a comparison of three cryptic binding site datasets: the PocketMiner test set, Cryptosite, and CryptoBench. Key statistics include average pocket RMSD, average number of binding residues per protein, and average number of observed residues per protein.
3. `stats.ipynb`  – Computes general statistics for the CryptoBench dataset, including the number of cryptic pockets, promiscuous pockets, and multi-chain pockets. It also calculates the number of apo-holo pairs remaining after each filtering phase,
4. `conservation` directory – Compares conservation score from the PDB API for between cryptic binding residues, non-cryptic binding residues, and non-binding residues,
5. `dataset-similarity` - Checks for similar structures between PocketMiner train set and CryptoBench test set.

### G-download-uniprot-sequences-and-create-annotations
This step prepares data for training and evaluating the benchmark method, organized into two subdirectories:
1. `download-uniprot-sequences` – Downloads UniProt sequences associated with the proteins in the CryptoBench dataset,
2. `create-annotations` – Extracts annotations from the CryptoBench dataset JSON file, generating the necessary annotation files for training and evaluation.

### Benchmark method training and evaluation
The training and evaluation of the benchmark method are handled by a framework stored in a [different repository](https://github.com/skrhakv/apolo/tree/cryptobench-v2), where you can find the complete workflow for these processes.

### H-prediction-evaluation
Although the evaluation of the benchmark method is managed in a separate repository, a distinct evaluation was required for PocketMiner. Therefore, the following notebooks are provided:
1. `roc-curve.ipynb` – Calculates the ROC and PRC curves to compare the performance of the benchmark method with PocketMiner,
2.  `skipped-structures-analysis.ipynb` – Analyzes any structures that were not evaluated by each method, identifying gaps in evaluation coverage.

### I-final
This folder contains supplementary materials, including example files, README documentation, and other data to be attached to the final dataset release.

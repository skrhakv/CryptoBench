import numpy as np
from typing import Dict
import csv
import statistics
import os
from sklearn import metrics
import pickle
import os

import biotite.structure.io.pdb as pdb
import biotite.structure as struc
from biotite.sequence import ProteinSequence

TEST_SUBSET_PATH = '/home/vit/Projects/cryptobench/data/G-download-uniprot-sequences-and-create-annotations/ahoj-v2/create-annotations/test.csv'
OUTPUT_PATH = '/home/vit/Projects/cryptobench/data/H-prediction-evaluation/ahoj-v2/pocketminer'
STATS_OUTPUT_PATH = f'{OUTPUT_PATH}/pocketminer-stats'
PDB_FILES = f'{OUTPUT_PATH}/pdb-cryptobench-v2-test-subset'
POCKETMINER_PREDICTIONS_PATH = f'{OUTPUT_PATH}/pocketminer-predictions'

class Results:
    def __init__(self, actual_values, predictions, predictions_for_auc, protein_code):
        self.actual_values = actual_values
        self.predictions = predictions
        self.cf = self.get_conf_matrix()
        self.predictions_for_auc = predictions_for_auc
        self.auc = self.get_auc()
        self.accuracy = self.get_accuracy()
        self.mcc = self.get_mcc()
        self.f1 = self.get_f1()


    def get_conf_matrix(self):
        return metrics.confusion_matrix(self.predictions, self.actual_values)

    def get_auc(self):
        return metrics.roc_auc_score(self.actual_values, self.predictions_for_auc)

    def get_accuracy(self):
        return metrics.accuracy_score(self.predictions, self.actual_values)

    def get_mcc(self):
        return metrics.matthews_corrcoef(self.predictions, self.actual_values)

    def get_f1(self):
        return metrics.f1_score(self.predictions, self.actual_values)

    def get_TPR(self):
        # FP = self.cf[1][0]
        FN = self.cf[0][1]
        TP = self.cf[1][1]
        # TN = self.cf[0][0]
        # Sensitivity, hit rate, recall, or true positive rate
        return TP / (TP + FN)

    def get_FPR(self):
        FP = self.cf[1][0]
        # FN = self.cf[0][1]
        # TP = self.cf[1][1]
        TN = self.cf[0][0]
        # Fall out or false positive rate
        return FP / (FP + TN)

def read_predictions() -> Dict[str, Results]:
    result = {}

    # using auth_seq_id annotations for the residue numbers and auth_chain_id for the chains
    with open(TEST_SUBSET_PATH, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            pdb_id = row[0].lower()
            chain_id = row[1].upper()
            annotations = row[3]

            id = f'{pdb_id}{chain_id}'

            # PocketMiner yields error for some structures -> they might not get predicted
            if not os.path.exists(f'{POCKETMINER_PREDICTIONS_PATH}/{pdb_id}{chain_id}.npy'):
                print(id)
                continue
            
            predictions = np.load(
                f'{POCKETMINER_PREDICTIONS_PATH}/{pdb_id}{chain_id}.npy')[0]
            pred_for_auc = np.copy(predictions)
            predictions[predictions > 0.5] = 1
            predictions[predictions <= 0.5] = 0

            annotations_set = set([i.split('_')[1] for i in annotations.split(' ')])

            cif_file = pdb.PDBFile.read(
                f'{PDB_FILES}/{pdb_id}{chain_id}.pdb')
            whole_structure = pdb.get_structure(cif_file, model=1)
            protein = whole_structure[struc.filter_amino_acids(
                whole_structure)]
            auth_ca_atoms = protein[(protein.atom_name == "CA") &
                                    (protein.element == "C")]
            y = [0] * predictions.shape[0]
            ix = 0

            for auth_ca in auth_ca_atoms:
                if auth_ca.chain_id != chain_id:
                    continue

                res_id = str(auth_ca).split()[1]
                if res_id in annotations_set:
                    y[ix - 1] = 1
                ix += 1

            # check that all binding residues are observed 
            assert sum(y) == len(annotations_set)
            result[f'{pdb_id}{chain_id}'] = Results(
                y, predictions, pred_for_auc, f'{pdb_id}{chain_id}')

    return result


def save_to_csv(result_list: Dict[str, Results], path):
    """Put everything into files"""
    with open(path, 'w', newline='') as file:
        writer = csv.writer(file)

        field = ["code", "length", "binding_residues",
                 "FPR", "TPR", "ACC", "MCC", "F1", "AUC"]

        writer.writerow(field)

        counter = 0
        avg_len = []
        avg_binding_res = []
        avg_fpr = []
        avg_tpr = []
        avg_acc = []
        avg_mcc = []
        avg_f1 = []
        avg_auc = []
        #
        # averaged OVER PROTEINS, NOT RESIDUES !!!
        #
        for protein_id, protein in result_list.items():
            counter += 1
            avg_len.append(len(protein.actual_values))
            avg_binding_res.append(sum(protein.actual_values))
            avg_fpr.append(protein.get_FPR())
            avg_tpr.append(protein.get_TPR())
            avg_acc.append(protein.accuracy)
            avg_mcc.append(protein.mcc)
            avg_f1.append(protein.f1)

            avg_auc.append(protein.auc)
            writer.writerow([protein_id, len(protein.actual_values), sum(protein.actual_values), protein.get_FPR(
            ), protein.get_TPR(), protein.accuracy, protein.mcc, protein.f1, protein.auc])

        writer.writerow(["average", sum(avg_len) / counter, sum(avg_binding_res) / counter, sum(avg_fpr) / counter, sum(
            avg_tpr) / counter, sum(avg_acc) / counter, sum(avg_mcc) / counter, sum(avg_f1) / counter, sum(avg_auc) / counter])
        writer.writerow(["standard deviation", statistics.stdev(avg_len), statistics.stdev(avg_binding_res), statistics.stdev(avg_fpr), statistics.stdev(
            avg_tpr), statistics.stdev(avg_acc), statistics.stdev(avg_mcc), statistics.stdev(avg_f1), statistics.stdev(avg_auc)])


test_set = read_predictions()
save_to_csv(test_set, f'{STATS_OUTPUT_PATH}/pocketminer-test-subset.csv')
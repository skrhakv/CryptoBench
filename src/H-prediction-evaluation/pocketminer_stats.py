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
VALIDATION_SUBSET_PATH = '/home/vit/Projects/cryptobench/data/G-download-uniprot-sequences-and-create-annotations/ahoj-v2/create-annotations/train-fold-0.csv'

OUTPUT_PATH = '/home/vit/Projects/cryptobench/data/H-prediction-evaluation/ahoj-v2/pocketminer'
STATS_OUTPUT_PATH = f'{OUTPUT_PATH}/pocketminer-stats'

TEST_PDB_FILES = f'{OUTPUT_PATH}/pdb-cryptobench-v2-test-subset'
VALIDATION_PDB_FILES = f'{OUTPUT_PATH}/pdb-cryptobench-v2-fold-0'

POCKETMINER_PREDICTIONS_PATH = f'{OUTPUT_PATH}/pocketminer-predictions'
POCKETMINER_VALIDATION_PREDICTIONS_PATH = f'{OUTPUT_PATH}/pocketminer-validation-predictions'
THRESHOLD = 0.75

mapping = {'Aba': 'A', 'Ace': 'X', 'Acr': 'X', 'Ala': 'A', 'Aly': 'K', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cas': 'C',
           'Ccs': 'C', 'Cme': 'C', 'Csd': 'C', 'Cso': 'C', 'Csx': 'C', 'Cys': 'C', 'Dal': 'A', 'Dbb': 'T', 'Dbu': 'T',
           'Dha': 'S', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'Glz': 'G', 'His': 'H', 'Hse': 'S', 'Ile': 'I', 'Leu': 'L',
           'Llp': 'K', 'Lys': 'K', 'Men': 'N', 'Met': 'M', 'Mly': 'K', 'Mse': 'M', 'Nh2': 'X', 'Nle': 'L', 'Ocs': 'C',
           'Pca': 'E', 'Phe': 'F', 'Pro': 'P', 'Ptr': 'Y', 'Sep': 'S', 'Ser': 'S', 'Thr': 'T', 'Tih': 'A', 'Tpo': 'T',
           'Trp': 'W', 'Tyr': 'Y', 'Unk': 'X', 'Val': 'V', 'Ycm': 'C', 'Sec': 'U', 'Pyl': 'O', 'Mhs': 'H', 'Snm': 'S',
           'Mis': 'S'}
def three_to_one(three_letter_code):
    return mapping[three_letter_code[0].upper() + three_letter_code[1:].lower()]


class Results:
    def __init__(self, actual_values, predictions, predictions_for_auc, protein_code):
        self.actual_values = actual_values
        self.predictions = predictions
        self.predictions_for_auc = predictions_for_auc
        self.cf = self.get_conf_matrix()
        self.auc = self.get_auc()
        self.accuracy = self.get_accuracy()
        self.mcc = self.get_mcc()
        self.f1 = self.get_f1()


    def get_conf_matrix(self):
        return metrics.confusion_matrix(self.predictions, self.actual_values)

    def get_auc(self):
        fpr, tpr, _ = metrics.roc_curve(self.actual_values, self.predictions_for_auc)
        roc_auc = metrics.auc(fpr, tpr)
        return roc_auc    

    def get_auprc(self):
        precision, recall, _ = metrics.precision_recall_curve(self.actual_values, self.predictions_for_auc)
        auc_precision_recall = metrics.auc(recall, precision)
        return auc_precision_recall

    def get_accuracy(self):
        return metrics.accuracy_score(self.predictions, self.actual_values)

    def get_mcc(self):
        return metrics.matthews_corrcoef(self.predictions, self.actual_values)

    def get_f1(self):
        return metrics.f1_score(self.predictions, self.actual_values, average='weighted')

    def satisfies_treshold(self, treshold_percent):
        return (metrics.confusion_matrix(self.predictions, self.actual_values)[1][1] / np.sum(self.actual_values) >= treshold_percent)

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

skipped_structures = []

def read_predictions(predictions_path, subset_path, pdb_files, threshold=0.5) -> Dict[str, Results]:
    result = {}

    # using auth_seq_id annotations for the residue numbers and auth_chain_id for the chains
    with open(subset_path, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            pdb_id = row[0].lower()
            chain_id = row[1].upper()
            annotations = row[3]

            id = f'{pdb_id}{chain_id}'

            if not os.path.exists(f'{predictions_path}/{pdb_id}{chain_id}.npy'):
                skipped_structures.append(id)
                print(id)
                continue
            if pdb_id == '7qni': continue

            predictions = np.load(
                f'{predictions_path}/{pdb_id}{chain_id}.npy')[0]
            pred_for_auc = np.load(
                f'{predictions_path}/{pdb_id}{chain_id}.npy')[0]
            tmp = predictions
            predictions[tmp > threshold] = 1
            predictions[tmp <= threshold] = 0

            annotations_set = set([i.split('_')[1] for i in annotations.split(' ')])

            cif_file = pdb.PDBFile.read(
                f'{pdb_files}/{pdb_id}{chain_id}.pdb')
            whole_structure = pdb.get_structure(cif_file, model=1)
            protein = whole_structure[struc.filter_amino_acids(
                whole_structure)]
            auth_ca_atoms = protein[(protein.atom_name == "CA") &
                                    (protein.element == "C")]
            y = [0] * predictions.shape[0]
            ix = 0

            if id == '7nbc':
                chain_id = 'C'
            for auth_ca in auth_ca_atoms:
                if auth_ca.chain_id != chain_id:
                    continue

                res_id = str(auth_ca).split()[1]
                if res_id in annotations_set:
                    y[ix] = 1
                ix += 1
            if id == '7nbc':
                chain_id = 'CCC'
            # check that all binding residues are observed 
            assert sum(y) == len(annotations_set), pdb_id
            result[f'{pdb_id}{chain_id}'] = Results(
                y, predictions, pred_for_auc, f'{pdb_id}{chain_id}')

            # print(y - predictions)
    return result


def save_to_csv(result_list: Dict[str, Results], path):
    with open(path, 'w', newline='') as file:
        writer = csv.writer(file)

        field = ["code", "length", "binding_residues",
                 "FPR", "TPR", "ACC", "MCC", "F1", "AUC", "AUPRC"]

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
        avg_auprc = []

        overall_predictions = None
        overall_predictions_for_auc = None
        overall_actual_values = None

        for protein_id, protein in result_list.items():

            if overall_predictions is None:
                overall_predictions = protein.predictions
                overall_predictions_for_auc = protein.predictions_for_auc
                overall_actual_values = protein.actual_values
            else:
                overall_predictions = np.concatenate((overall_predictions, protein.predictions), axis=0)
                overall_predictions_for_auc = np.concatenate((overall_predictions_for_auc, protein.predictions_for_auc), axis=0)
                overall_actual_values = np.concatenate((overall_actual_values, protein.actual_values), axis=0)

            counter += 1
            avg_len.append(len(protein.actual_values))
            avg_binding_res.append(sum(protein.actual_values))
            avg_fpr.append(protein.get_FPR())
            avg_tpr.append(protein.get_TPR())
            avg_acc.append(protein.accuracy)
            avg_mcc.append(protein.mcc)
            avg_f1.append(protein.f1)
            avg_auprc.append(protein.get_auprc())

            avg_auc.append(protein.auc)
            writer.writerow([protein_id, len(protein.actual_values), sum(protein.actual_values), protein.get_FPR(
            ), protein.get_TPR(), protein.accuracy, protein.mcc, protein.f1, protein.auc, protein.get_auprc()])
        overall = Results(overall_actual_values, overall_predictions, overall_predictions_for_auc, 'overall')

        writer.writerow(["average", sum(avg_len) / counter, sum(avg_binding_res) / counter, sum(avg_fpr) / counter, sum(
            avg_tpr) / counter, sum(avg_acc) / counter, sum(avg_mcc) / counter, sum(avg_f1) / counter, sum(avg_auc) / counter, sum(avg_auprc) / counter])
        writer.writerow(["standard deviation", statistics.stdev(avg_len), statistics.stdev(avg_binding_res), statistics.stdev(avg_fpr), statistics.stdev(
            avg_tpr), statistics.stdev(avg_acc), statistics.stdev(avg_mcc), statistics.stdev(avg_f1), statistics.stdev(avg_auc), statistics.stdev(avg_auprc)])
        writer.writerow(["overall", overall_predictions.shape[0], sum(overall.actual_values), overall.get_FPR(), overall.get_TPR(),overall.accuracy, overall.mcc, overall.f1, overall.auc, overall.get_auprc() ])

def to_labels(pos_probs, threshold):
    return (pos_probs >= threshold).astype('int')

def get_optimal_decision_threshold():
    result_list = read_predictions(POCKETMINER_VALIDATION_PREDICTIONS_PATH, VALIDATION_SUBSET_PATH, VALIDATION_PDB_FILES)
    for protein in result_list.values():
        overall_probs = None
        overall_actual = None

        if overall_probs is None:
            overall_probs = protein.predictions_for_auc
            overall_actual = protein.actual_values
        else:
            overall_probs = np.concatenate((overall_probs, protein.predictions_for_auc), axis=0)
            overall_actual = np.concatenate((overall_actual, protein.actual_values), axis=0)
        
    thresholds = np.arange(0, 1, 0.05)
        
    fpr, tpr, thresholds = metrics.roc_curve(overall_actual, overall_probs)
    print(thresholds)
    scores = [metrics.f1_score(overall_actual, to_labels(overall_probs, t), average='weighted') for t in thresholds]
    # get the best threshold
    J = tpr + scores - fpr
    ix = np.argmax(J)
    print('Threshold=%.3f' % (thresholds[ix]))
    return thresholds[ix]

if __name__ == '__main__':
    threshold = get_optimal_decision_threshold()
    print(f'threshold is {threshold} ..')
    test_set = read_predictions(POCKETMINER_PREDICTIONS_PATH, TEST_SUBSET_PATH, TEST_PDB_FILES, threshold)   

    # with open('results.pickle', 'wb') as handle:
    #     pickle.dump(test_set, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # with open('results.pickle', 'rb') as handle:
    #     test_set = pickle.load(handle)
    #
    # os.makedirs('stats', exist_ok=True)
    os.makedirs('stats', exist_ok=True)
    save_to_csv(test_set, f'{STATS_OUTPUT_PATH}/pocketminer-test-subset.csv')   

    # save skipped structures
    with open(f'{OUTPUT_PATH}/skipped-structures.txt', 'w') as f:
        for i in skipped_structures:
            f.write(f'{i}\n')
import numpy as np
from sklearn import metrics  

class Protein:
    def __init__(self, id, predictions, actual_values, threshold=0.5):
        self.id: str = id
        self.actual_values: np.ndarray = actual_values

        self.predictions = np.copy(predictions[:, 1])
        self.predictions[self.predictions > threshold] = 1
        self.predictions[self.predictions <= threshold] = 0

        self.cf = self.get_conf_matrix()
        self.auc = self.get_auc(predictions)
        self.auprc = self.get_auprc(predictions)
        self.accuracy = self.get_accuracy()
        self.mcc = self.get_mcc()
        self.f1 = self.get_f1()
    
    def get_conf_matrix(self):
        return metrics.confusion_matrix(self.predictions, self.actual_values)

    def get_auc(self, pred):
        fpr, tpr, _ = metrics.roc_curve(self.actual_values, pred[:,1])
        roc_auc = metrics.auc(fpr, tpr)
        return roc_auc    

    def get_auprc(self, pred):
        precision, recall, _ = metrics.precision_recall_curve(self.actual_values, pred[:,1])
        auc_precision_recall = metrics.auc(recall, precision)
        return auc_precision_recall

    def get_accuracy(self):
        return metrics.accuracy_score(self.predictions, self.actual_values)
    
    def get_mcc(self):
        return metrics.matthews_corrcoef(self.predictions, self.actual_values)
    
    def get_f1(self):
        return metrics.f1_score(self.predictions, self.actual_values)
    
    def satisfies_treshold(self, treshold_percent):
        return (metrics.confusion_matrix(self.predictions, self.actual_values)[1][1] / np.sum(self.actual_values) >= treshold_percent)
    
    def get_TPR(self):
        # FP = self.cf[1][0]
        FN = self.cf[0][1]
        TP = self.cf[1][1]
        # TN = self.cf[0][0]
        # Sensitivity, hit rate, recall, or true positive rate
        assert FN + TP == np.sum(self.actual_values)
        return TP/(TP+FN)
    
    def get_FPR(self):
        FP = self.cf[1][0]
        # FN = self.cf[0][1]
        # TP = self.cf[1][1]
        TN = self.cf[0][0]
        # Fall out or false positive rate
        assert FP + TN == len(self.actual_values) - np.sum(self.actual_values)
        return FP/(FP+TN)

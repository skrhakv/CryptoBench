import numpy as np
from Protein import Protein
from tensorflow import keras
import tensorflow_addons as tfa

# example script to predict cryptic pockets for the '7w19A' apo structure (from CryptoBench test set)

MODEL_PATH = '../benchmark/best_trained'
STRUCTURE_ID = '7w19A'

# 0.75 decision threshold was used in the CryptoBench paper
DECISION_THRESHOLD = 0.75


def load_model():
    print("Loading CryptoBench model ...")
    return keras.models.load_model(MODEL_PATH,
                                   custom_objects={
                                       'MatthewsCorrelationCoefficient': tfa.metrics.MatthewsCorrelationCoefficient(num_classes=2)},
                                   compile=False)


def predict(X, model):
    print("Making prediction ...")
    return model.predict(X)


def load_data():
    print("Loading data - embeddings and annotations ...")
    embeddings = np.load(f'data/{STRUCTURE_ID}.npy')

    #
    # CAUTION: If you compare the annotation in 'data/annotation.txt' and in the original
    # '../cryptobench-dataset/dataset.json', the annotation for 7w19A is shifted by 2.
    # There are two reasons for this:
    # 1. The first residue is not observed, therefore it was excluded during extraction of the embedding
    # 2. The author labeling at 7w19A starts at one. However, annotation in 'annotation.txt' is zero-based.
    # 
    with open('data/annotation.txt', 'r') as f:
        annotations = f.read().split(';')[3].split(' ')

    # the format of each annotation is as follows: 
    # 'A_N22' denotes a single binding residue, which belongs to the 'A' chain,
    # 'N' denotes that the residue is Asparagine, and the corresponding embedding
    # can be found at index 22
    annotations = [int(i.split('_')[1][1:]) for i in annotations]
    y = [0] * embeddings.shape[0]
    for ix in annotations:
        y[ix] = 1

    return embeddings, y


def print_evaluation(evaluation):
    print(
        f'\n\n\nEvaluation for {evaluation.id} with decision threshold = {DECISION_THRESHOLD}:\n')
    print(f'AUC: {evaluation.auc}')
    print(f'AUPRC: {evaluation.auprc}')
    print(f'ACC: {evaluation.accuracy}')
    print(f'TPR: {evaluation.get_TPR()}')
    print(f'FPR: {evaluation.get_FPR()}')
    print(f'MCC: {evaluation.mcc}')
    print(f'F1: {evaluation.f1}')


def evaluate(prediction, actual):
    evaluation = Protein(STRUCTURE_ID, prediction, actual,
                         threshold=DECISION_THRESHOLD)
    print_evaluation(evaluation)


def main():
    model = load_model()
    embeddings, annotations = load_data()
    predictions = predict(embeddings, model)
    evaluate(predictions, annotations)


if __name__ == '__main__':
    main()

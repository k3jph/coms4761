import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn import metrics


def main():
    print("Loading in data...")
    data_pos = pd.read_csv('../data/datapos.csv', header=None)
    data_neg = pd.read_csv('../data/datapos.csv', header=None)

    pos_train, pos_test = train_test_split(data_pos,test_size=0.3)
    neg_train, neg_test = train_test_split(data_neg,test_size=0.3)

    train = np.concatenate((pos_train, neg_train))
    train_labels = np.concatenate((np.ones(pos_train.shape[0]), np.zeros(neg_train.shape[0])))

    test = np.concatenate((pos_test, neg_test))
    test_labels = np.concatenate((np.ones(pos_test.shape[0]), np.zeros(neg_test.shape[0])))

    print("Training model...")
    log_reg = LogisticRegression(penalty='l1', solver='liblinear')
    log_reg.fit(train, train_labels)

    label_preds = log_reg.predict_proba(test)[::,0]

    auc = metrics.roc_auc_score(test_labels, label_preds)
    print(f"Area Under Curve (AUC): {auc}")


if __name__ == '__main__':
    main()
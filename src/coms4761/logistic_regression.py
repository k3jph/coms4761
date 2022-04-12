import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.metrics import plot_precision_recall_curve, plot_roc_curve, roc_auc_score


def main():
    print("Loading in data...")
    data_pos = pd.read_csv('../data/datapos.csv', header=None)
    data_neg = pd.read_csv('../data/datapos.csv', header=None)

    pos_train, pos_test = train_test_split(data_pos, test_size=0.1, shuffle=True)
    neg_train, neg_test = train_test_split(data_neg, test_size=0.1, shuffle=True)

    train = np.concatenate((pos_train, neg_train))
    train_labels = np.concatenate((np.ones(pos_train.shape[0]), np.zeros(neg_train.shape[0])))

    test = np.concatenate((pos_test, neg_test))
    test_labels = np.concatenate((np.ones(pos_test.shape[0]), np.zeros(neg_test.shape[0])))

    cv = KFold(n_splits=7, random_state=1, shuffle=True)

    print("Training LR model...")
    log_reg = LogisticRegression(penalty="l1", solver="liblinear")

    scores = cross_val_score(log_reg, train, train_labels, cv=cv)
    print(f"Cross Validation Scores: {scores}")

    log_reg.fit(train, train_labels)

    label_preds = log_reg.predict(test)

    display = plot_precision_recall_curve(log_reg, test, test_labels, name="Linear")
    _ = display.ax_.set_title("2-class Precision-Recall curve")
    plt.show()

    auc_display = plot_roc_curve(log_reg, test, test_labels)
    plt.show()

    auc = roc_auc_score(test_labels, label_preds)
    print(f"Area Under Curve (AUC): {auc}") 

    print("Training RF model...")
    rf = RandomForestClassifier()
    rf.fit(train, train_labels)

    label_preds = rf.predict_proba(test)
    predictions = rf.predict(test)

    display = plot_precision_recall_curve(rf, test, test_labels, name="RF")
    _ = display.ax_.set_title("2-class Precision-Recall curve")
    plt.show()

    auc_display = plot_roc_curve(rf, test, test_labels)
    plt.show()

    auc = roc_auc_score(predictions, test_labels)
    print(f"Area Under Curve (AUC): {auc}") 


if __name__ == '__main__':
    main()
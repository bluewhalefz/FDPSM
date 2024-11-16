import argparse

import joblib
import pandas as pd
import numpy as np
import os
import time
from openfe import OpenFE, transform
from sklearn.metrics import classification_report, confusion_matrix, recall_score, precision_score, f1_score, \
    matthews_corrcoef, roc_auc_score, precision_recall_curve, auc, roc_curve, accuracy_score
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def predict(train_file, test_file, autoFE_features, selection, model_file):
    train = pd.read_csv(train_file, low_memory=False)
    test = pd.read_csv(test_file, low_memory = False)
    fea_list = test.columns.tolist()
    fea_list = [ele.replace("-", "_") for ele in fea_list]
    test.columns = fea_list
    test_info = test.iloc[:, :10]

    features = joblib.load(autoFE_features)
    

    X_train = train.iloc[:, 10:].astype('float64')
    X_test = test.iloc[:, 10:].astype('float64')
    print(X_test.shape)

    print('---' + time.asctime(time.localtime(time.time())) + '--- transforming dataset\n')
    _, X_test_tr = transform(X_train, X_test, features, n_jobs=30)
    feature_list_final = pd.read_csv(selection)["feature"].tolist()
    print(feature_list_final)
    X_test_filtered = X_test_tr[feature_list_final].astype('float64')
    print(X_test_filtered)

    print('---' + time.asctime(time.localtime(time.time())) + '--- predicting\n')
    model = joblib.load(model_file)
    test_pred = model.predict(X_test_filtered)
    return pd.concat([test_info, X_test_filtered, pd.DataFrame(test_pred, columns=['FDPSM'])], axis=1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_file", type=str, help="输入文件")
    parser.add_argument("--test_file", type=str, help="输入文件")
    parser.add_argument("--openfe_features")
    parser.add_argument("--selectedFea")
    parser.add_argument("--model_file", type=str)
    parser.add_argument("--outfile", type=str, help="输出文件")

    args = parser.parse_args()

    # test_file = "/data1/jinfangfang/Project/ODDSM/data/test816_hg38_98fea_processed.csv"
    # autoFE_features = "/data1/jinfangfang/Project/ODDSM/result/openFE_train8502_hg38_98fea_processed.features"
    # selection = "/data1/jinfangfang/Project/ODDSM/result/selection_train8502_hg38_98fea_processed.csv"
    # model_file = "/data1/jinfangfang/Project/ODDSM/result/train8502_hg38_98fea_processed.model"
    # filename = "test816_98fea"
    df_pred = predict(args.train_file, args.test_file, args.openfe_features, args.selectedFea, args.model_file)
    df_pred.to_csv(args.outfile, index=False)
    name = "FDPSM"
    df_pred["y_pred_classes"] = np.where(df_pred[name] > 0.65, 1, 0)
    y_label = df_pred["label"].tolist()
    y_pred = df_pred[name].tolist()
    y_pred_classes = df_pred["y_pred_classes"].tolist()
    classification_metrics = classification_report(y_label, y_pred_classes)
    print(classification_metrics)

    cm = confusion_matrix(y_label, y_pred_classes)
    tn, fp, fn, tp = cm.ravel()

    # SEN
    sen = recall_score(y_label, y_pred_classes)

    # SPE
    spe = tn / (tn + fp)

    # PRE
    pre = precision_score(y_label, y_pred_classes)

    # F1
    f1 = f1_score(y_label, y_pred_classes)

    #  MCC
    mcc = matthews_corrcoef(y_label, y_pred_classes)

    print("SEN = %.4f" % sen)
    print("SPE = %.4f" % spe)
    print("PRE = %.4f" % pre)
    print("F1 = %.4f" % f1)
    print("MCC = %.4f" % mcc)

    auc_score = roc_auc_score(y_label, y_pred)

    precision, recall, thresholds = precision_recall_curve(y_label, y_pred)
    aupr = auc(recall, precision)
    print("AUC = %.4f" % auc_score, "    AUPR =%.4f" % aupr)

    fpr, tpr, thresholds = roc_curve(y_label, y_pred)
    auc_test = auc(fpr, tpr)

    acc = accuracy_score(y_label, y_pred_classes)
    print("ACC = %.4f" % acc)




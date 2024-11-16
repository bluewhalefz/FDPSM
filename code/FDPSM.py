import argparse
import os.path
import joblib
import numpy as np
import pandas as pd
import time
import re
import itertools as it

from lightgbm import early_stopping
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from openfe import OpenFE, transform, tree_to_formula
import lightgbm as lgb
import warnings


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

selection_params = {
    'boosting': 'dart',
    'objective': 'binary',
    'metric': {'binary_logloss', 'auc'},  # 二进制对数损失
    'drop_rate':0.1,
    'drop_seed':1,
    'max_drop':5,
    'num_leaves': 31,
    'max_depth':5,     #ver1：max_depth : 6   AUC为0.73是max_depth为5，min_data_in_leaf为51，num_leaves为30
    'max_bin': 255,
    'min_data_in_leaf': 61,    #ver：min_data_in_leaf = 101
    'learning_rate': 0.01,    # ver1、ver2:'learning_rate': 0.01
    'feature_fraction': 1.0,
    'bagging_fraction': 1.0,
    'bagging_freq': 45,
    'lambda_l1': 0.001,
    'lambda_l2': 0.4,  # 越小l2正则程度越高
    'min_split_gain': 0.0,
    'verbose': 5,
    'is_unbalance': False
}

threshold = 0.5
beta_score = 0.5

callbacks = [early_stopping(stopping_rounds=10)]
callbacks1 = [early_stopping(stopping_rounds=5)]
callbacks2 = [early_stopping(stopping_rounds=50)]

def load_data(df_train, y_train):
    train_x, train_y = df_train, y_train
    X, val_X, y, val_y = train_test_split(
        train_x,
        train_y,
        test_size=0.2,
        random_state=1,
        stratify=train_y
    )
    lgb_train = lgb.Dataset(X, y)
    lgb_eval = lgb.Dataset(val_X, val_y, reference=lgb_train)
    return lgb_train, lgb_eval

def train_model(df_train, y_train):
    lgb_train, lgb_eval = load_data(df_train, y_train)
    gbm = lgb.train(selection_params,
                    lgb_train,
                    num_boost_round=40,
                    valid_sets=[lgb_eval],
                    callbacks=callbacks1
                    )
    return gbm

def autofe(data, openfe_features, selectedFea_file, importance_file, feaName_file):
    print('---' + time.asctime(time.localtime(time.time())) + '--- autoFE\n')
    # train_info记录前四列的信息，Chrom、Pos、Ref Alt CLASS
    train_info = data.iloc[:, :10]
    X_train = data.iloc[:, 10:].astype('float64')
    _, X_test = train_test_split(X_train, test_size=0.1)
    Y_train = np.array(data['label']).reshape(len(data['label']), )

    ofe = OpenFE()
    # generate new features
    features = ofe.fit(data=X_train, label=Y_train, n_jobs=30)
    # 执行特征变换：transform()方法可能会对原始数据集中的特征进行各种转换操作，如缩放、归一化、标准化、离散化等。这些变换可以是基于之前
    # fit()方法中学到的统计信息，也可以是一些预定义的转换规则。
    X_train_tr, X_test_tr = transform(X_train, X_test, features, n_jobs=30)
    print("**********************finished***********************")
    # joblib.dump(features, os.path.join(root, f'result/openFE_{filename}.features'))
    joblib.dump(features, openfe_features)
    #
    X_train_tr.index = list(range(X_train_tr.shape[0]))
    
    with open(feaNama_file, 'w') as f_write:
        f_write.write('name\n')
        for feature in features:
            f_write.write(tree_to_formula(feature) + '\n')

    
    feature_num = X_train_tr.shape[1]
    print("feature_num:", feature_num)
    train_round = 1

    while feature_num > 200:
        gbm = train_model(X_train_tr, Y_train)
        feature_imp = pd.DataFrame({'Value': gbm.feature_importance(), 'Feature': X_train_tr.columns})

        if train_round >= 40:
            break
        train_round += 1
        feature_sum = feature_imp.Value.sum()
        drop_list = feature_imp[feature_imp.Value / feature_sum <= 0].Feature.tolist()

        df = feature_imp
        df["importance"] = feature_imp.Value / feature_sum
        df = df[df["importance"] > 0]
        df.to_csv(importance_file, index=False)

        X_train_tr.drop(drop_list, axis=1, inplace=True)
        feature_num = X_train_tr.shape[1]
    print("feature_num: ", feature_num)
    df = pd.DataFrame(X_train_tr.columns, columns=['feature'])
    df.to_csv(selectedFea_file, index=False)
    return pd.concat([train_info, X_train_tr], axis=1)
    


def step_training(data, model_file):
    X_train = data.iloc[:, 10:]
    Y_train = data['label']
    X_train = X_train.astype("float64")
    lgb_train = lgb.Dataset(X_train, Y_train, free_raw_data=False)

    # set init params excluding CV params
    params = {
        'boosting_type': 'gbdt',
        'objective': 'binary',
        'metric': 'auc',
        'nthread': 4,
        'learning_rate': 0.1
    }
    max_auc = float('0')
    best_params = {}

    # Imporve accuracy
    for num_leaves in range(5, 100, 5):
        for max_depth in range(3, 8, 1):
            params['num_leaves'] = num_leaves
            params['max_depth'] = max_depth

            cv_results = lgb.cv(
                params,
                lgb_train,
                seed=1,
                nfold=5,
                metrics=['auc'],
                callbacks=callbacks,
                eval_train_metric=True
            )

            mean_auc = pd.Series(cv_results['valid auc-mean']).max()

            if mean_auc >= max_auc:
                max_auc = mean_auc
                best_params['num_leaves'] = num_leaves
                best_params['max_depth'] = max_depth
    if 'num_leaves' and 'max_depth' in best_params.keys():
        params['num_leaves'] = best_params['num_leaves']
        params['max_depth'] = best_params['max_depth']

    # Avoid over-fitting
    for max_bin in range(5, 256, 10):
        for min_data_in_leaf in range(1, 102, 10):
            params['max_bin'] = max_bin
            params['min_data_in_leaf'] = min_data_in_leaf

            cv_results = lgb.cv(
                params,
                lgb_train,
                seed=1,
                nfold=5,
                metrics=['auc'],
                callbacks=callbacks,
                eval_train_metric=True
            )

            mean_auc = pd.Series(cv_results['valid auc-mean']).max()

            if mean_auc >= max_auc:
                max_auc = mean_auc
                best_params['max_bin'] = max_bin
                best_params['min_data_in_leaf'] = min_data_in_leaf
    if 'max_bin' and 'min_data_in_leaf' in best_params.keys():
        params['min_data_in_leaf'] = best_params['min_data_in_leaf']
        params['max_bin'] = best_params['max_bin']

    for feature_fraction in [0.6, 0.7, 0.8, 0.9, 1.0]:
        for bagging_fraction in [0.6, 0.7, 0.8, 0.9, 1.0]:
            for bagging_freq in range(0, 50, 5):
                params['feature_fraction'] = feature_fraction
                params['bagging_fraction'] = bagging_fraction
                params['bagging_freq'] = bagging_freq

                cv_results = lgb.cv(
                    params,
                    lgb_train,
                    seed=1,
                    nfold=5,
                    metrics=['auc'],
                    callbacks=callbacks,
                    eval_train_metric=True
                )

                mean_auc = pd.Series(cv_results['valid auc-mean']).max()
                boost_rounds = pd.Series(cv_results['valid auc-mean']).idxmax()

                if mean_auc >= max_auc:
                    max_auc = mean_auc
                    best_params['feature_fraction'] = feature_fraction
                    best_params['bagging_fraction'] = bagging_fraction
                    best_params['bagging_freq'] = bagging_freq

    if 'feature_fraction' and 'bagging_fraction' and 'bagging_freq' in best_params.keys():
        params['feature_fraction'] = best_params['feature_fraction']
        params['bagging_fraction'] = best_params['bagging_fraction']
        params['bagging_freq'] = best_params['bagging_freq']

    for lambda_l1 in [1e-5, 1e-3, 1e-1, 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
        for lambda_l2 in [1e-5, 1e-3, 1e-1, 0.0, 0.1, 0.4, 0.6, 0.7, 0.9, 1.0]:
            params['lambda_l1'] = lambda_l1
            params['lambda_l2'] = lambda_l2
            cv_results = lgb.cv(
                params,
                lgb_train,
                seed=1,
                nfold=5,
                metrics=['auc'],
                callbacks=callbacks,
                # verbose_eval=True
                eval_train_metric=True
            )

            mean_auc = pd.Series(cv_results['valid auc-mean']).max()
            boost_rounds = pd.Series(cv_results['valid auc-mean']).idxmax()

            if mean_auc >= max_auc:
                max_auc = mean_auc
                best_params['lambda_l1'] = lambda_l1
                best_params['lambda_l2'] = lambda_l2
    if 'lambda_l1' and 'lambda_l2' in best_params.keys():
        params['lambda_l1'] = best_params['lambda_l1']
        params['lambda_l2'] = best_params['lambda_l2']

    for min_split_gain in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        params['min_split_gain'] = min_split_gain

        cv_results = lgb.cv(
            params,
            lgb_train,
            seed=1,
            nfold=5,
            metrics=['auc'],
            callbacks=callbacks,
            # verbose_eval=True
            eval_train_metric=True
        )

        mean_auc = pd.Series(cv_results['valid auc-mean']).max()
        boost_rounds = pd.Series(cv_results['valid auc-mean']).idxmax()

        if mean_auc >= max_auc:
            max_auc = mean_auc

            best_params['min_split_gain'] = min_split_gain
    if 'min_split_gain' in best_params.keys():
        params['min_split_gain'] = best_params['min_split_gain']

    print("**********params*********")
    print(params)
    print("*************************")
    lgb_train, lgb_eval = load_data(X_train, Y_train)
    params = {
        'boosting_type': 'gbdt',
        'objective': 'binary',
        'metric': {'binary_logloss', 'auc'},  # 二进制对数损失
        'num_leaves': 40,
        'max_depth': 6,
        'max_bin': 255,
        'min_data_in_leaf': 101,
        'learning_rate': 0.01,    
        'feature_fraction': 1.0,
        'bagging_fraction': 1.0,
        'bagging_freq': 45,
        'lambda_l1': 0.001,
        'lambda_l2': 0.4,  # 越小l2正则程度越高
        'min_split_gain': 0.0,
        'verbose': 5,
        'is_unbalance': False
    }
    for key in best_params.keys():
        if key == 'max_depth':
            params[key] = best_params[key]
        elif key == 'max_leaves':
            params[key] = best_params[key]
        else:
            params[key] = best_params[key]

    gbm = lgb.train(params,
                    lgb_train,
                    num_boost_round=1000,
                    valid_sets=lgb_eval,
                    callbacks=callbacks2)
    joblib.dump(gbm, os.path.join(root, model_file))




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_file", type=str, help="inputfile: training set")
    parser.add_argument("--selectedFea")
    parser.add_argument("--openfe_features")
    parser.add_argument("--model_file", type=str, help="输出文件")
    parser.add_argument("--importance_file", type=str, help="输出文件")
    parser.add_argument("--feaName_file", type=str, help="输出文件")

    args = parser.parse_args()
    data = pd.read_csv(args.train_file)
    data = autofe(data, args.openfe_features, args.selectedFea, args.importance_file, args.feaName_file)
    #step_training(data, args.model_file)
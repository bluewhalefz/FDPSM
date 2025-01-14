import argparse
import os.path
import joblib
import numpy as np
import pandas as pd
import time
import re
import itertools as it

from lightgbm import early_stopping
from sklearn.model_selection import train_test_split, RandomizedSearchCV,StratifiedKFold
from sklearn.metrics import roc_auc_score, make_scorer
from openfe import OpenFE, transform, tree_to_formula
import lightgbm as lgb
import warnings
import os
from scipy.stats import uniform, randint, loguniform
import random



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
    'max_depth':5,     #ver1：max_depth : 6   
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
        random_state=2024,
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
   
    joblib.dump(features, openfe_features)
    
    X_train_tr.index = list(range(X_train_tr.shape[0]))
    
    with open(feaName_file, 'w') as f_write:
        f_write.write('name\n')
        for feature in features:
            f_write.write(tree_to_formula(feature) + '\n')

    
    feature_num = X_train_tr.shape[1]
    print("feature_num:", feature_num)
    train_round = 0

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
    


def step_training(data, model_file, FDPSM_importance):
    X = data.iloc[:, 10:]
    Y = data['label']
    X = X.astype("float64")
    X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size=0.2, random_state=2024)
     

    param_dist = {
    'num_leaves': range(5, 100, 5),
    'max_depth': range(3, 8, 1),
    'max_bin': range(5, 256, 10),
    'min_data_in_leaf': range(1, 102, 10),
    'learning_rate': [0.001, 0.004, 0.007, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4,0.5],
    #'learning_rate': np.logspace(-3, 0, 10),  # 在0.001到1之间对数均匀分布
    'feature_fraction': np.linspace(0.6, 1.0, 5),  # 在0.6到1之间线性均匀分布
    'bagging_fraction': np.linspace(0.6, 1.0, 5),
    'bagging_freq': range(0, 50, 5),
    'lambda_l1': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    'lambda_l2': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    #'lambda_l2': np.logspace(-5, 0, 6),
    'min_split_gain': np.linspace(0.0, 1.0, 11)
    }

    fixed_params = {
        'boosting_type': 'gbdt',
        'objective': 'binary',
        'metric': 'auc',
        'nthread': 4,
        'verbose': -1  # 设置为-1以减少输出，或在调试时设置为更高的值
    }
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=2024)
    
    
    # 初始化LGBMClassifier
    lgb_clf = lgb.LGBMClassifier(**fixed_params)
    # 使用RandomizedSearchCV进行随机搜索和嵌套交叉验证
    random_search = RandomizedSearchCV(
        estimator=lgb_clf,
        param_distributions=param_dist,
        n_iter=5000,  # 尝试的参数组合数量
        scoring='roc_auc',  # 使用AUC作为评分指标
        cv=cv,  # 使用分层交叉验证
        n_jobs=-1,  # 使用所有可用的CPU核心
        verbose=1,  # 显示进度
        random_state=2024,  # 设置随机种子以确保结果可重复
    )
    
    random_search.fit(X, Y)
    
    # 输出最优参数
    print("**********Best Params*********")
    print(random_search.best_params_)
    print("*************************")
    print(f"Best AUC: {random_search.best_score_:.4f}")
    
    lgb_train, lgb_eval = load_data(X, Y) 
    
    # 使用最佳参数训练最终模型
    best_params = random_search.best_params_
    final_params = {**fixed_params, **best_params}
    final_params['metric'] = ['binary_logloss', 'auc']  # 更新为所需的评估指标
    gbm = lgb.train(
        final_params,
        lgb_train,
        num_boost_round=1000,  # 
        valid_sets=lgb_eval,
        callbacks=callbacks2
        #verbose=True  # 每100次迭代输出一次结果（根据需要调整）
    )
    
    feature_imp = pd.DataFrame({'Value': gbm.feature_importance(), 'Feature': X.columns})
    feature_sum = feature_imp.Value.sum()
    df = feature_imp
    df["importance"] = feature_imp.Value / feature_sum
    df.to_csv(FDPSM_importance, index=False)
     
    # 保存模型
    joblib.dump(gbm, os.path.join(root, model_file))  # 确保root变量已定义


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_file", type=str, help="输入文件")
    parser.add_argument("--selectedFea")
    parser.add_argument("--openfe_features")
    parser.add_argument("--model_file", type=str, help="输出文件")
    parser.add_argument("--importance_file", type=str, help="输出文件")
    parser.add_argument("--feaName_file", type=str, help="输出文件")
    parser.add_argument("--FDPSM_importance", type=str)

    args = parser.parse_args()
    data = pd.read_csv(args.train_file)
    #data = autofe(data, args.openfe_features, args.selectedFea, args.importance_file, args.feaName_file)
    print(data.shape)
    step_training(data, args.model_file, args.FDPSM_importance)
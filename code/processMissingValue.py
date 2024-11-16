import argparse
import pandas as pd
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge


# 使用基于贝叶斯模型的迭代式方法对缺失值进行填充
def processMissing(train_file, test_file, out_train, out_test):
    df_train = pd.read_csv(train_file)
    listFea_train = df_train.columns.tolist()
    colnames_train = listFea_train[10:]


    df_test = pd.read_csv(test_file)
    listFea_test = df_test.columns.tolist()
    new_listFea_test = listFea_test[:10] + colnames_train
    df_test = df_test.reindex(columns=new_listFea_test)

    # 实例化一个BatesianRidge模型作为迭代填充器
    estimator = BayesianRidge()

    # 实例化IterativeImputer
    imputer = IterativeImputer(estimator=estimator, random_state=2024, max_iter=10)

    # 使用贝叶斯PCA填充缺失值
    train_filled = imputer.fit_transform(df_train.iloc[:, 10:])
    test_filled = imputer.transform(df_test.iloc[:, 10:])

    # 将填充后的数据转换为df
    
    train_filled_df = pd.DataFrame(train_filled, columns=colnames_train)
    print("特征的数量一共为（不包括chr、POS、ref、ALT）", train_filled_df.shape[1])
    train_filled_df = pd.concat([df_train.iloc[:, :10], train_filled_df], axis=1)

    test_filled_df = pd.DataFrame(test_filled, columns=colnames_train)
    print("特征的数量一共为（不包括chr、POS、ref、ALT）", test_filled_df.shape[1])
    test_filled_df = pd.concat([df_test.iloc[:, :10], test_filled_df], axis=1)
    
    
    train_filled_df.to_csv(out_train, index=False)
    test_filled_df.to_csv(out_test, index=False)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_file", type=str, help="训练集")
    parser.add_argument("--test_file", type=str, help="测试集")
    parser.add_argument("--outfile_train", type=str, help="经过缺失值填充之后的测试集")
    parser.add_argument("--outfile_test", type=str, help="经过缺失值填充之后的测试集")

    args = parser.parse_args()
    processMissing(args.train_file, args.test_file, args.outfile_train, args.outfile_test)




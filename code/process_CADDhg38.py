import argparse

import pandas as pd
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge


def process(input, file, outfile):
    df = pd.read_csv(file)
    df["Chrom"] = df["Chrom"].astype(str)
    df["Pos"] = df["Pos"].astype(str)
    print(len(df))
    # 首先将列“ConsDetail”中的注释包含”synonymous“的行全部提取，保存到df_synonymous
    df_synonymous = df[df["ConsDetail"].str.contains("synonymous")]
    print(len(df_synonymous))
    #然后将原文件删除df_synonymous，保存到df_resiual，并将df_resiual进行去重（按照"#Chrom", "Pos", "Ref", "Alt"这四列进行去重）
    df_resiual = pd.concat([df, df_synonymous], axis=0)
    df_resiual.drop_duplicates(subset=None, keep=False, inplace=True)
    # print(len(df_resiual))
    #
    df_resiual = df_resiual.sort_values("PHRED", ascending=False)
    df_resiual.drop_duplicates(subset=["Chrom", "Pos", "Ref", "Alt"], keep="first", inplace=True)
    # print(len(df_resiual))
    #最后将df_synonymous和df_resiual进行行拼接，之后进行去重操作（按照"#Chrom", "Pos", "Ref", "Alt"这四列进行去重）
    df_synonymous = pd.concat([df_synonymous, df_resiual], axis=0)
    df_synonymous = df_synonymous.sort_values("PHRED", ascending=False)
    df_synonymous.drop_duplicates(subset=["Chrom", "Pos", "Ref", "Alt"], keep="first", inplace=True)

    print(len(df_synonymous))
    df_all = pd.read_csv(input, header=None, sep="\t", usecols=[0, 1, 3, 4])
    df_all.columns=["Chrom", "Pos", "Ref", "Alt"]
    df_all["Chrom"] = df_all["Chrom"].astype(str)
    df_all["Pos"] = df_all["Pos"].astype(str)
    df_merged = pd.merge(df_all, df_synonymous, how="left", on=["Chrom", "Pos", "Ref", "Alt"])
    print(df_merged.shape[1])
    df_merged.drop(columns=["ConsDetail", "PHRED"], inplace=True)
    print(df_merged.shape[1])
    df_merged.to_csv(outfile, sep=",", index=False)

# 并将CADD与已注释好的特征拼接起来
def mergefea(caddfile, annotated, allfea_file):
    df_cadd = pd.read_csv(caddfile)
    df_annotated = pd.read_csv(annotated)

    df_cadd["CHROM_hg19"] = df_cadd["CHROM_hg19"].astype(str)
    df_cadd["POS_hg19"] = df_cadd["POS_hg19"].astype(str)

    df_cadd["CHROM_hg38"] = df_cadd["CHROM_hg38"].astype(str)
    df_cadd["POS_hg38"] = df_cadd["POS_hg38"].astype(str)

    df_annotated["Chr"] = df_annotated["Chr"].astype(str)
    df_annotated["Start"] = df_annotated["Start"].astype(str)


    # 注意：有annotated的注释信息是基于hg38的
    df_annotated["CHROM_hg38"] = df_annotated["Chr"].str.split("r", expand=True)[1]
    df_annotated.rename(columns={"Start":"POS_hg38", "Ref":"REF", "Alt":"ALT"}, inplace=True)
    df_merged = pd.merge(df_cadd, df_annotated, on=["CHROM_hg38", "POS_hg38", "REF", "ALT"], how="left")
    df_merged.to_csv(allfea_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfile", help="原始所有突变输入文件的路径")
    parser.add_argument("--cadd_file", help="具有cadd注释的突变文件")
    parser.add_argument("--outfile", help="输出文件的路径")

    args = parser.parse_args()

    process(args.inputfile, args.cadd_file, args.outfile)


    # annotated_train = dir + "Annotation/train8502_hg38_out.csv"
    # annotated_test = dir + "Annotation/test816_hg38_out.csv"
    # allfea_train = dir + "Annotation/train8502_feature_hg38.csv"
    # allfea_test = dir + "Annotation/test816_feature_hg38.csv"
    # mergefea(out_train, annotated_train, allfea_train)
    # mergefea(out_test, annotated_test, allfea_test)


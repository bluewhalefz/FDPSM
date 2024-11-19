import argparse

import pandas as pd
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge


def process(inputfile, source_file, outfile):
    df = pd.read_csv(inputfile)
    df.rename(columns = {"Chrom":"CHROM_hg38", "Pos":"POS_hg38"}, inplace=True)
    df["CHROM_hg38"] = "chr" + df["CHROM_hg38"]
    df["POS_hg38"] = df["POS_hg38"].astype(str)
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
    df_resiual.drop_duplicates(subset=["CHROM_hg38", "POS_hg38", "Ref", "Alt"], keep="first", inplace=True)
    # print(len(df_resiual))
    #最后将df_synonymous和df_resiual进行行拼接，之后进行去重操作（按照"#Chrom", "Pos", "Ref", "Alt"这四列进行去重）
    df_synonymous = pd.concat([df_synonymous, df_resiual], axis=0)
    df_synonymous = df_synonymous.sort_values("PHRED", ascending=False)
    df_synonymous.drop_duplicates(subset=["CHROM_hg38", "POS_hg38", "Ref", "Alt"], keep="first", inplace=True)

    print(len(df_synonymous))
    df_all = pd.read_csv(source_file)
    df_all.rename(columns = {"Chr":"CHROM_hg38", "Start":"POS_hg38"}, inplace=True)
    df_all["CHROM_hg38"] = df_all["CHROM_hg38"].astype(str)
    df_all["POS_hg38"] = df_all["POS_hg38"].astype(str)
    df_merged = pd.merge(df_all, df_synonymous, how="left", on=["CHROM_hg38", "POS_hg38", "Ref", "Alt"])
    print(df_merged.shape[1])
    df_merged.drop(columns=["ConsDetail", "PHRED"], inplace=True)
    print(df_merged.shape[1])
    colnames = df_merged.columns.tolist()
    colnames = [name.replace("-", "_") for name in colnames]
    df_merged.columns = colnames
    fea_list = ["CHROM_hg38", "POS_hg38", "Ref", "Alt",  'Dist2Mutation', 'Freq100bp', 'Rare100bp', 'Sngl100bp', 'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp', 'Rare10000bp', 'Sngl10000bp', 'priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP', 'verPhyloP', 'bStatistic', 'ZooPriPhyloP', 'ZooVerPhyloP', 'GC', 'CpG', 'minDistTSS', 'minDistTSE', 'GerpN', 'GerpS', 'RegSeq0', 'RegSeq1', 'RegSeq2', 'RegSeq3', 'RegSeq4', 'RegSeq5', 'RegSeq6', 'RegSeq7', 'cHmm_E1', 'cHmm_E2', 'cHmm_E3', 'cHmm_E4', 'cHmm_E5', 'cHmm_E6', 'cHmm_E7', 'cHmm_E8', 'cHmm_E9', 'cHmm_E10', 'cHmm_E11', 'cHmm_E12', 'cHmm_E13', 'cHmm_E14', 'cHmm_E15', 'cHmm_E16', 'cHmm_E17', 'cHmm_E18', 'cHmm_E19', 'cHmm_E20', 'cHmm_E21', 'cHmm_E22', 'cHmm_E23', 'cHmm_E24', 'cHmm_E25', 'EncodeH3K4me1_sum', 'EncodeH3K4me1_max', 'EncodeH3K4me2_sum', 'EncodeH3K4me2_max', 'EncodeH3K4me3_sum', 'EncodeH3K4me3_max', 'EncodeH3K9ac_sum', 'EncodeH3K9ac_max', 'EncodeH3K9me3_sum', 'EncodeH3K9me3_max', 'EncodeH3K27ac_sum', 'EncodeH3K27ac_max', 'EncodeH3K27me3_sum', 'EncodeH3K27me3_max', 'EncodeH3K36me3_sum', 'EncodeH3K36me3_max', 'EncodeH3K79me2_sum', 'EncodeH3K79me2_max', 'EncodeH4K20me1_sum', 'EncodeH4K20me1_max', 'EncodeH2AFZ_sum', 'EncodeH2AFZ_max', 'EncodeDNase_sum', 'EncodeDNase_max']
    df_final = df_merged[fea_list].reindex(columns = fea_list)
    df_final.to_csv(outfile, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfile")
    parser.add_argument("--source_file")
    parser.add_argument("--outfile", help="输出文件的路径")

    args = parser.parse_args()

    process(args.inputfile, args.source_file, args.outfile)




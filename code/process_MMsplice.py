import argparse

import pandas as pd

def process(inputfile, test_file, outfile):
    df = pd.read_csv(inputfile, sep="\t")
    df[["num@CHROM_hg38", "POS_hg38", "Ref", "Alt"]] = df["key"].str.split("_", expand=True)
    df[["num", "CHROM_hg38"]] = df["num@CHROM_hg38"].str.split("@", expand=True)
    df.drop(columns=["num", "num@CHROM_hg38", "key"], inplace=True)

    df["CHROM_hg38"] = "chr" + df["CHROM_hg38"]
    df["POS_hg38"] = df["POS_hg38"].astype(str)
    df_test = pd.read_csv(test_file)
    df_test["CHROM_hg38"] = df_test["CHROM_hg38"].astype(str)
    df_test["POS_hg38"] = df_test["POS_hg38"].astype(str)

    df_match = pd.merge(df_test, df, how="left", on=["CHROM_hg38", "POS_hg38", "Ref", "Alt"])
    colnames = df_match.columns.tolist()
    colnames = [name.replace("-", "_") for name in colnames]
    df_match.columns = colnames
    
    #df = pd.read_csv("/data1/jinfangfang/Project/ODDSM/data/train8502_hg38_111fea.csv")
    #fea_list = df.iloc[:, 10:].columns.tolist()
    fea_list = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'MMS_delta_logit_psi', 'MMS_ref_acceptorIntron', 'MMS_alt_acceptorIntron', 'MMS_ref_acceptor', 'MMS_alt_acceptor', 'MMS_ref_exon', 'MMS_alt_exon', 'MMS_ref_donor', 'MMS_alt_donor', 'MMS_ref_donorIntron', 'MMS_alt_donorIntron', 'MMS_pathogenicity', 'MMS_efficiency', 'Dist2Mutation', 'Freq100bp', 'Rare100bp', 'Sngl100bp', 'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp', 'Rare10000bp', 'Sngl10000bp', 'priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP', 'verPhyloP', 'bStatistic', 'ZooPriPhyloP', 'ZooVerPhyloP', 'GC', 'CpG', 'minDistTSS', 'minDistTSE', 'GerpN', 'GerpS', 'RegSeq0', 'RegSeq1', 'RegSeq2', 'RegSeq3', 'RegSeq4', 'RegSeq5', 'RegSeq6', 'RegSeq7', 'gdi', 'gdi_phred', 'cHmm_E1', 'cHmm_E2', 'cHmm_E3', 'cHmm_E4', 'cHmm_E5', 'cHmm_E6', 'cHmm_E7', 'cHmm_E8', 'cHmm_E9', 'cHmm_E10', 'cHmm_E11', 'cHmm_E12', 'cHmm_E13', 'cHmm_E14', 'cHmm_E15', 'cHmm_E16', 'cHmm_E17', 'cHmm_E18', 'cHmm_E19', 'cHmm_E20', 'cHmm_E21', 'cHmm_E22', 'cHmm_E23', 'cHmm_E24', 'cHmm_E25', 'EncodeH3K4me1_sum', 'EncodeH3K4me1_max', 'EncodeH3K4me2_sum', 'EncodeH3K4me2_max', 'EncodeH3K4me3_sum', 'EncodeH3K4me3_max', 'EncodeH3K9ac_sum', 'EncodeH3K9ac_max', 'EncodeH3K9me3_sum', 'EncodeH3K9me3_max', 'EncodeH3K27ac_sum', 'EncodeH3K27ac_max', 'EncodeH3K27me3_sum', 'EncodeH3K27me3_max', 'EncodeH3K36me3_sum', 'EncodeH3K36me3_max', 'EncodeH3K79me2_sum', 'EncodeH3K79me2_max', 'EncodeH4K20me1_sum', 'EncodeH4K20me1_max', 'EncodeH2AFZ_sum', 'EncodeH2AFZ_max', 'EncodeDNase_sum', 'EncodeDNase_max', 'Hmec', 'Hsmm', 'Huvec', 'K562', 'Nhek', 'Nhlf']
   
    
    
    df_fea = df_match[fea_list].reindex(columns = fea_list)
    df_match = pd.concat([df_match.iloc[:, :10], df_fea], axis = 1)
    df_match.to_csv(outfile, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfile", type=str, help="输入文件")
    parser.add_argument("--source", type=str, help="测试集文件")
    parser.add_argument("--outfile", type=str, help="输出文件")

    args = parser.parse_args()

    process(args.inputfile, args.source, args.outfile)


import argparse

import pandas as pd

def process(inputfile, source_file, outfile):
    df = pd.read_csv(inputfile, sep="\t")
    df[["num@CHROM_hg38", "POS_hg38", "Ref", "Alt"]] = df["key"].str.split("_", expand=True)
    df[["num", "CHROM_hg38"]] = df["num@CHROM_hg38"].str.split("@", expand=True)
    df.drop(columns=["num", "num@CHROM_hg38", "key"], inplace=True)

    df["CHROM_hg38"] = "chr" + df["CHROM_hg38"]
    df["POS_hg38"] = df["POS_hg38"].astype(str)
    df_source = pd.read_csv(source_file)
    df_source["CHROM_hg38"] = df_source["CHROM_hg38"].astype(str)
    df_source["POS_hg38"] = df_source["POS_hg38"].astype(str)

    df_match = pd.merge(df_source, df, how="left", on=["CHROM_hg38", "POS_hg38", "Ref", "Alt"])

    df_match.to_csv(outfile, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfile", type=str, help="输入文件")
    parser.add_argument("--source", type=str, help="原始文件")
    parser.add_argument("--outfile", type=str, help="输出文件")

    args = parser.parse_args()

    process(args.inputfile, args.source, args.outfile)


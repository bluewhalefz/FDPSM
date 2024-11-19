from datetime import datetime
import pandas as pd
import numpy as np
import vcf
import time
import re
import os
import argparse
import paramiko
import os

SPLICEAI_DATA_DIR = "./data/spliceai/"
ANNOVAR_DATA_DIR = "./data/annovar/"
ANNOVAR_DIR = "./annovar/"


spliceai_column_list = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']
root = "./"
epi_file = os.path.join(root, 'data/other_database/master38.chromhmm.bedg')
rvis_file = os.path.join(root, 'data/other_database/RVIS.txt')
gdi_file = os.path.join(root, 'data/other_database/GDI.txt')
lof_file = os.path.join(root, 'data/other_database/Lof.txt')


chr_list = ['chr' + str(chromosome) for chromosome in range(1, 23)] + ['chrX', 'chrY']
func_dict = {'nonsynonymous SNV': 'nonsynonymous SNV', 'synonymous SNV': 'synonymous SNV',
             'stopgain': 'stopgain', 'unknown': 'unknown', 'startloss': 'startloss',
             'stoploss': 'stoploss', 'frameshift substitution': 'frameshift',
             'frameshift insertion': 'frameshift', 'frameshift deletion': 'frameshift',
             'nonframeshift insertion': 'nonframeshift', 'nonframeshift deletion': 'nonframeshift',
             'nonframeshift substitution': 'nonframeshift',
             'NA': 'NA'}
type_dict = {'0': 0,
             'Active_Promoter': 1,
             'Heterochrom/lo': 2,
             'Insulator': 3,
             'Poised_Promoter': 4,
             'Repetitive/CNV': 5,
             'Repressed': 6,
             'Strong_Enhancer': 7,
             'Txn_Elongation': 8,
             'Txn_Transition': 9,
             'Weak_Enhancer': 10,
             'Weak_Promoter': 11,
             'Weak_Txn': 12,
             0: 0}


def prepare_input(input_file, filename, annovar_dir, spliceai_dir):
    """
    Prepare input file for ANNOVAR and SpliceAI.
    Input file should be tab or comma separated files containing five headers: Chr, Start, End, Ref, Alt
    :param annovar_dir: Dictionary of ANNOVAR input file. (e.g. /data/annovar/)
    :param spliceai_dir: Dictionary of SpliceAI input file. (e.g. /data/spliceai/)
    :return: None
    """
    print(input_file)
    data = pd.read_csv(input_file, low_memory=False)
    # print(data)
    data.Chr = [str(i).replace('chr', '') if 'chr' in str(i) else str(i) for i in data.Chr.tolist()]
    data = data[['Chr', 'Start', 'End', 'Ref', 'Alt']]
    # data = data[data.Chr.isin(chr_list)]
    data = data.fillna('-')
    data.insert(data.shape[-1], 'info', '.')
    print(data)
    data.to_csv(f'{os.path.join(annovar_dir, filename)}.avinput', sep='\t', index=False, header=None)
    data.to_csv(f'{os.path.join(spliceai_dir, filename)}.vcf', sep='\t', index=False, header=None)
    now = datetime.now()
    str_to_insert = f'''##fileformat=VCFv4.1
##fileDate={"{:02d}-{:02d}-{:02d}".format(now.year, now.month, now.day)}
##source=magpie
##reference=GRCh38
##contig=<ID=1,length=248956422>
##contig=<ID=2,length=242193529>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##ID=<Description="ClinVar Variation ID">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
'''
    temp_filename = f"{os.path.join(spliceai_dir, filename)}_temp.txt"
    spliceai_file = f'{os.path.join(spliceai_dir, filename)}.vcf'
    with open(temp_filename, 'w') as temp_file:
        temp_file.write(str_to_insert)
        with open(spliceai_file, 'r') as file:
            temp_file.write(file.read())
    os.remove(spliceai_file)
    os.rename(temp_filename, spliceai_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfile", type=str)
    parser.add_argument("--filename", type=str)
    args = parser.parse_args()
    print("Received inputfile:", args.inputfile)

    prepare_input(args.inputfile.strip(), args.filename, ANNOVAR_DATA_DIR, SPLICEAI_DATA_DIR)





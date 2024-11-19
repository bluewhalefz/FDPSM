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

def binarySearch(list1, target):
    left = 0
    right = len(list1) - 1

    while left <= right:
        mid = (right + left) // 2
        if mid + 1 < len(list1) and list1[mid] <= target < list1[mid + 1]:
            return mid
        elif target > list1[mid]:
            left = mid + 1
        else:  # target < list[mid]
            right = mid - 1
    return -1



def annotate_spliceai(df, vcf_file):
    print('---' + time.asctime(time.localtime(time.time())) + '--- annotating SpliceAI\n')
    vcf_data = vcf.Reader(filename = vcf_file)
    data_lists = []
    index = 0
    gene_list = df['Gene.refGene'].to_list()
    for record in vcf_data:
        temp_list = []
        if record.INFO.get('SpliceAI') is not None:
            spliceai_list = record.INFO.get('SpliceAI')[0].split(',')
            for spliceai in spliceai_list:
                temp_list = spliceai.split('|')
                if temp_list[1] == gene_list[index]:
                    temp_list = [float(i) if i != '.' else 0.0 for i in temp_list[2:]]
                else:
                    temp_list = [0.0] * 8
        else:
            temp_list = [0.0] * 8
        data_lists.append(temp_list)
        index += 1
    if df.shape[0] == len(data_lists):
        data_spliceai = pd.DataFrame(np.matrix(data_lists))
        data_spliceai.columns = spliceai_column_list
        df = df.join(data_spliceai, how='left')
        return df
    else:
        assert 'Cannot add features to data'


def annotate_chromHMM(data):  # after annotating spliceai
    print('---' + time.asctime(time.localtime(time.time())) + '--- annotating chromHMM.\n')

    cell_types = ['Gm12878', 'H1hesc', 'Hepg2', 'Hmec', 'Hsmm', 'Huvec', 'K562', 'Nhek', 'Nhlf']

    data.Chr = ['chr' + str(i) for i in data.Chr.tolist() if 'chr' not in str(i)]
    data = data[data.Chr.isin(chr_list)]
    data = data.sort_values(by=['Chr', 'Start', 'End', 'Ref', 'Alt'], ascending=True)

    # containing position
    chr_dict = {}
    # containing content
    content_dict = {}
    for c in chr_list:
        chr_dict[c] = []
        content_dict[c] = []
    for c in chr_list:
        count = 0
        with open(epi_file, 'r') as f_read:
            temp = f_read.readline()
            while True:
                line = f_read.readline()
                if line:
                    line = line.strip().split('\t')
                    if line[0] == str(c) and count == 0:
                        chr_dict[c].extend([int(line[1]), int(line[2])])
                        content_dict[c].append(line[3:])
                    elif line[0] == str(c) and count != 0:
                        chr_dict[c].append(int(line[2]))
                        content_dict[c].append(line[3:])
                    count += 1
                else:
                    break

    epi_list = []
    for c in chr_list:
        for record in data[data.Chr == c].iterrows():
            binary_start_index = binarySearch(chr_dict[c], int(record[1]['Start']))
            end_index = binarySearch(chr_dict[c], int(record[1]['End']))
            if binary_start_index == -1 and end_index == -1:
                epi_list.append([0] * 9)
            elif binary_start_index != -1 and end_index == -1:
                epi_list.append(content_dict[c][binary_start_index])
            elif binary_start_index == -1 and end_index != -1:
                epi_list.append(content_dict[c][end_index])
            elif binary_start_index == end_index:
                epi_list.append(content_dict[c][end_index])
            elif binary_start_index != end_index:
                temp = []
                for index in range(9):
                    if content_dict[c][binary_start_index][index] == content_dict[c][end_index][index]:
                        temp.append(content_dict[c][binary_start_index][index])
                    else:
                        temp.append(0)
                epi_list.append(temp)

    s = set()
    for epi in epi_list:
        s.add(len(epi))
#    print(len(s))
#    print(len(epi_list) == data.shape[0])

    epi_dict = {}
    for epi in cell_types:
        epi_dict[epi] = []
    for index in range(len(epi_list)):
        count = 0
        for epi in cell_types:
            epi_dict[epi].append(epi_list[index][count])
            count += 1

    for col in cell_types:
        data.insert(data.shape[-1], col, epi_dict[col])

    for col in cell_types:
        data[col] = [type_dict[cell_type] for cell_type in data[col].to_list()]

    data = data.sort_index(ascending = True)
    data.index = list(range(data.shape[0]))

    return data

def annotate_functional_effect(data):
    print('---' + time.asctime(time.localtime(time.time())) + '--- ' + 'annotating functional effects\n')

    # data['ExonicFunc.refGene'].fillna('unknown', inplace=True)
    # func_list = data['ExonicFunc.refGene'].tolist()
    # func_list = [func_dict[i] for i in func_list]
    # data.drop('ExonicFunc.refGene', axis=1, inplace=True)
    # data.insert(data.shape[-1], 'func', func_list)
    #
    gene_list = data['Gene.refGene'].unique().tolist()
    sum_gdi = {}
    sum_gdi_phred = {}
    all_disease = {}
    all_Mendelian = {}
    mendelian_AD = {}
    mendelian_AR = {}
    all_PID = {}
    pid_AD = {}
    pid_AR = {}
    count = 0
    with open(gdi_file, 'r') as read_GDI:
        while True:
            line = read_GDI.readline()
            if line:
                list_text = line.split()
                if str(list_text[0]) in gene_list:
                    sum_gdi[list_text[0]] = float(list_text[1])
                    sum_gdi_phred[list_text[0]] = float(list_text[2])
                    all_disease[list_text[0]] = list_text[3]
                    all_Mendelian[list_text[0]] = list_text[4]
                    mendelian_AD[list_text[0]] = list_text[5]
                    mendelian_AR[list_text[0]] = list_text[6]
                    all_PID[list_text[0]] = list_text[7]
                    pid_AD[list_text[0]] = list_text[8]
                    pid_AR[list_text[0]] = list_text[9]
                    count += 1
                if count == len(gene_list):
                    break
            else:
                break

    for i in gene_list:
        if i not in sum_gdi.keys():
            sum_gdi[i] = np.nan
            sum_gdi_phred[i] = np.nan
            all_disease[i] = np.nan
            all_Mendelian[i] = np.nan
            mendelian_AD[i] = np.nan
            mendelian_AR[i] = np.nan
            all_PID[i] = np.nan
            pid_AD[i] = np.nan
            pid_AR[i] = np.nan

    count = 0
    rvis1 = {}
    rvis2 = {}
    with open(rvis_file, 'r') as read_RVIS:
        while True:
            line = read_RVIS.readline()
            if line:
                list_text = line.split()
                if list_text[0] in gene_list or list_text[1] in gene_list or \
                        list_text[2] in gene_list or list_text[3] in gene_list or \
                        list_text[4] in gene_list:
                    rvis1[list_text[0]] = float(list_text[5])
                    rvis2[list_text[0]] = float(list_text[6])
                    count += 1
                if count == len(gene_list):
                    break
            else:
                break

    for i in gene_list:
        if i not in rvis1.keys():
            rvis1[i] = np.nan
            rvis2[i] = np.nan

    count = 0
    lof_score = {}
    with open(lof_file, 'r') as read_lof:
        while True:
            line = read_lof.readline()
            if line:
                list_text = line.split()
                if str(list_text[0]) in gene_list:
                    lof_score[str(list_text[0])] = float(list_text[1])
                    count += 1
                if count == len(gene_list):
                    break
            else:
                break

    for i in gene_list:
        if i not in lof_score.keys():
            lof_score[i] = np.nan

    sum_gdi_list = []
    sum_gdi_phred_list = []
    rvis1_list = []
    rvis2_list = []
    lof_score_list = []
    for record in data['Gene.refGene']:
        sum_gdi_list.append(sum_gdi[record])
        sum_gdi_phred_list.append(sum_gdi_phred[record])
        rvis1_list.append(rvis1[record])
        rvis2_list.append(rvis2[record])
        lof_score_list.append(lof_score[record])

    lists = [sum_gdi_list, sum_gdi_phred_list, rvis1_list, rvis2_list, lof_score_list]
    column_name_list = ['gdi', 'gdi_phred', 'rvis1', 'rvis2', 'lof_score']
    column_name_dict = {}
    for col_index in range(len(column_name_list)):
        column_name_dict[col_index] = column_name_list[col_index]
    data_gene = pd.DataFrame(lists).T
    data_gene.rename(columns=column_name_dict, inplace=True)
    data = data.join(data_gene, how='left')
    return data

def select(data):
    feature_list = ['Chr', 'Start', 'Ref', 'Alt', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL','Hmec', 'Hsmm',
                    'Huvec', 'K562', 'Nhek', 'Nhlf','gdi', 'gdi_phred']
    data = data[feature_list]
    data.rename(columns = {"Chr":"CHROM_hg38", "Start":"POS_hg38"}, inplace=True)
    return data

def annotate(input_file, spliceai_output):
    data = pd.read_csv(input_file, low_memory = False)
    data = annotate_spliceai(data, spliceai_output)
    data = annotate_chromHMM(data)
    data = annotate_functional_effect(data)
    data = select(data)
    return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfile", type=str)
    parser.add_argument("--filename", type=str)
    parser.add_argument("--outfile", type=str)
    args = parser.parse_args()
    
    
    out_spliceai = SPLICEAI_DATA_DIR + args.filename + "_out.vcf"

    data = annotate(args.inputfile, out_spliceai)
    data.to_csv(args.outfile, index=False)




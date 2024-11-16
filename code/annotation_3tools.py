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

SPLICEAI_DATA_DIR = "/data1/jinfangfang/Project/MAGPIE/data/output/spliceai/"
ANNOVAR_DATA_DIR = "/data1/jinfangfang/Project/MAGPIE/data/output/annovar/"
ANNOVAR_DIR = "/data1/jinfangfang/Project/MAGPIE/annovar/"

class MyLinux(object):
    # 通过IP, 用户名，密码，超时时间初始化一个远程Linux主机
    def __init__(self, ip, username, password, timeout=30):
        self.ip = ip
        self.username = username
        self.password = password
        self.timeout = timeout
        # transport和chanel
        self.t = ''
        self.chan = ''
        # 链接失败的重试次数
        self.try_times = 3
        self.port = 22
        self.ssh_client = paramiko.SSHClient()

    # 调用该方法连接远程主机
    def connect(self):
        while True:
            # 连接过程中可能会抛出异常，比如网络不通、链接超时
            try:
                self.t = paramiko.Transport(sock=(self.ip, 22))
                self.t.connect(username=self.username, password=self.password)
                self.chan = self.t.open_session()
                self.chan.settimeout(self.timeout)
                self.chan.get_pty()
                self.chan.invoke_shell()
                # 如果没有抛出异常说明连接成功，直接返回
                print('连接%s成功' % self.ip)
                # 接收到的网络数据解码为str
                return
            except Exception as e1:
                if self.try_times != 0:
                    print('连接%s失败，进行重试' % self.ip)
                    self.try_times -= 1
                else:
                    print('重试3次失败，结束程序')
                    exit(1)

    # 断开连接
    def close(self):
        # 如果要运行服务，不能关闭连接
        self.chan.close()
        self.t.close()

    # 发送要执行的命令
    def send(self, cmd):
        cmd += '\r'
        # 发送要执行的命令
        self.chan.send(cmd)
        # 回显很长的命令可能执行较久，通过循环分批次取回回显,执行成功返回true,失败返回false
        result = ''
        while True:
            time.sleep(100)   #time.sleep(0.5)
            ret = self.chan.recv(65535)
            result += ret.decode('utf-8')
            if len(ret) < 65535:
                break
        print(result.replace(result.split("\n")[-1], "").strip("\n"))


spliceai_column_list = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']
root = "/data1/jinfangfang/Project/MAGPIE/"
epi_file = os.path.join(root, 'data/annotation_database/master38.chromhmm.bedg')
rvis_file = os.path.join(root, 'data/annotation_database/RVIS.txt')
gdi_file = os.path.join(root, 'data/annotation_database/GDI.txt')
lof_file = os.path.join(root, 'data/annotation_database/Lof.txt')


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

def prepare_input(input_file, annovar_dir, spliceai_dir):
    """
    Prepare input file for ANNOVAR and SpliceAI.
    :param input_file: Original input filename of MAGPIE. Input file should be tab or comma separated files containing five headers: Chr, Start, End, Ref, Alt
    :param annovar_dir: Dictionary of ANNOVAR input file. (e.g. /data/output/annovar/)
    :param spliceai_dir: Dictionary of SpliceAI input file. (e.g. /data/output/spliceai/)
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

def executeCommand(filename):
    host = MyLinux('172.17.4.67', 'jinfangfang', 'jinfangfang')  # 传入Ip，用户名，密码
    host.connect()

    annovar_path = ANNOVAR_DIR + "table_annovar.pl"
    data = ANNOVAR_DATA_DIR + filename + ".avinput"
    data_path = ANNOVAR_DATA_DIR + "humandb/"
    outfile = ANNOVAR_DATA_DIR + filename
    host.send("nohup perl {0} {1} {2} -buildver hg38 -out {3} -remove "
              "-protocol refGene,phastConsElements100way,gnomad30_genome,dbnsfp33a,dbnsfp42a "
              "-operation g,r,f,f,f -csvout &".format(annovar_path, data, data_path, outfile))
    time.sleep(2000)

    host.send("conda activate spliceai")

    spliceai_path = SPLICEAI_DATA_DIR + filename + ".vcf"
    out_spliceai_path = SPLICEAI_DATA_DIR + filename + "_out.vcf"
    spliceai_hg38 = SPLICEAI_DATA_DIR + "hg38.fa"
    host.send("CUDA_VISIBLE_DEVICES=2,1 nohup spliceai -I {0} -O {1} -R {2} -A grch38 &".format(spliceai_path, out_spliceai_path, spliceai_hg38))
    time.sleep(150000)
    host.send("conda deactivate")


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
    parser.add_argument("--inputfile", type=str, default=None)
    parser.add_argument("--outfile", type=str, default=None)
    args = parser.parse_args()

    filename = args.inputfile.split("/")[-1].split(".")[0]
    out_spliceai = SPLICEAI_DATA_DIR + filename + "_out.vcf"

    # 准备输入文件，将输入文件进行处理
    prepare_input(args.inputfile, ANNOVAR_DATA_DIR, SPLICEAI_DATA_DIR)
    executeCommand(filename)

    inputfile = ANNOVAR_DATA_DIR + filename + ".hg38_multianno.csv"
    data = annotate(inputfile, out_spliceai)
    data.to_csv(args.outfile, index=False)




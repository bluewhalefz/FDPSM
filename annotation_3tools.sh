#!/bin/bash
ANNOVAR_DIR="./annovar/"
ANNOVAR_DATA_DIR="./data/output/annovar/"
SPLICEAI_DATA_DIR="./data/output/spliceai/"
DATA_FILE="$1"
OUTPUT_FILE="$2"
FILE_NAME="$3"
#source activate magpie
#python code/prepare_input.py --inputfile "${DATA_FILE} --filename "${FILE_NAME}""
#echo "table_annovar.pl path is ${ANNOVAR_DIR}table_annovar.pl"
#echo "path:${FILE_NAME}"
#perl ${ANNOVAR_DIR}table_annovar.pl ${ANNOVAR_DATA_DIR}"${FILE_NAME}".avinput ${ANNOVAR_DATA_DIR}humandb/ -buildver hg38 -out ${ANNOVAR_DATA_DIR}"${FILE_NAME}" -remove -protocol refGene,phastConsElements100way,gnomad30_genome,dbnsfp33a,dbnsfp42a -operation g,r,f,f,f -csvout
#source deactivate
source activate spliceai
CUDA_VISIBLE_DEVICES=2,1 spliceai -I ${SPLICEAI_DATA_DIR}"${FILE_NAME}".vcf -O ${SPLICEAI_DATA_DIR}"${FILE_NAME}"_out.vcf -R ${SPLICEAI_DATA_DIR}hg38.fa -A grch38
conda deactivate
source activate magpie
python code/annotate.py --inputfile ${ANNOVAR_DATA_DIR}"${FILE_NAME}".hg38_multianno.csv --filename "${FILE_NAME}" --outfile "${OUTPUT_FILE}"

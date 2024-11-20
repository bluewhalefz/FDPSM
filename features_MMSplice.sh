#!/bin/bash

export hg19="data/Homo_sapiens/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
export hg38="data/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
export hg19_gtf="data/gtf/gencode.v34lift37.annotation.gtf"
export hg38_gtf="data/gtf/gencode.v34.annotation.gtf"


if [ $# -lt 3 ]; then
    echo "usage: $0 vcf buildver outdir"
    exit 1
fi
## input
input_vcf="$1"
export buildver="$2"
 if [ $buildver != "hg19"  ] && [ $buildver != "hg38" ]; then
#if [ $buildver != "hg19"  ]; then
    echo - "Unknown build version: $buildver (hg19 required)"
    exit 1
fi
export outdir="$3"

SRC="`dirname $0`/code"

bn=`basename $input_vcf`
bn=${bn%.gz}
export bn=${bn%.vcf}

export WORKING_DIR=$(dirname `realpath $0`)


## configuration of variables
if [ $buildver = "hg19" ]; then
    export ref_fasta="$hg19"
    export gtf="$hg19_gtf"
elif [ $buildver = "hg38" ]; then
    export ref_fasta="$hg38"
    export gtf="$hg38_gtf"
fi

## ======================= main ============================= 
## 1. format
vcf=${outdir}/${bn}.std.vcf
$SRC/format_vcf.py $input_vcf > $vcf
$SRC/vcf_add_header.py $vcf -r $ref_fasta --overwrite
avinput=${vcf%.vcf}.avinput
export source_file="$4"
export outfile="$5"

## mmsplice (hg19/hg38, dependent on reference genome and annotation)
if [ -e ${avinput}.mmsplice.tsv ]; then
    echo -e "- skip mmsplice"
else
    echo -e "- run mmsplice => ${avinput}.mmsplice.tsv"
    python $SRC/mmsplice_predict.py -vcf $vcf --gtf $gtf --fasta $ref_fasta -o ${avinput}.mmsplice.tsv  -bs  256 > ${avinput}.mmsplice.tsv.log
fi
$SRC/prepare_mmsplice.py -vcf $vcf -m ${avinput}.mmsplice.tsv 2> ${avinput}.features-mmsplice.tsv.log > ${avinput}.features-mmsplice.tsv || exit "ERROR in mmsplice module!"
python $SRC/process_MMsplice.py  --inputfile ${avinput}.features-mmsplice.tsv --source $source_file --outfile $outfile
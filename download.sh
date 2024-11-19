#fasta file
wget -O data/spliceai/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip data/spliceai/hg38.fa.gz

#chromHMM 15 core state
tar -zxvf data/annotation_database/master38.chromhmm.bedg.tar.gz

# AnnoVar annotation file
wget -O data/annovar/humandb/hg38_phastConsElements100way.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_phastConsElements100way.txt.gz
wget -O data/annovar/humandb/hg38_refGene.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_refGene.txt.gz
wget -O data/annovar/humandb/hg38_refGeneMrna.fa.gz http://www.openbioinformatics.org/annovar/download/hg38_refGeneMrna.fa.gz
wget -O data/annovar/humandb/hg38_refGeneVersion.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_refGeneVersion.txt.gz
wget -O data/annovar/humandb/hg38_dbnsfp33a.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp33a.txt.gz
wget -O data/annovar/humandb/hg38_dbnsfp33a.txt.idx.gz http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp33a.txt.idx.gz
wget -O data/annovar/humandb/hg38_dbnsfp41a.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp41a.txt.gz
wget -O data/annovar/humandb/hg38_dbnsfp41a.txt.idx.gz http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp41a.txt.idx.gz
wget -O data/annovar/humandb/hg38_dbnsfp42a.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp42a.txt.gz
wget -O data/annovar/humandb/hg38_gnomad30_genome.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_gnomad30_genome.txt.gz
wget -O data/annovar/humandb/hg38_gnomad30_genome.txt.idx.gz http://www.openbioinformatics.org/annovar/download/hg38_gnomad30_genome.txt.idx.gz
wget -O data/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -O data/Homo_sapiens/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz https://ftp.ensembl.org/pub/grch37/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
wget -O data/gtf/gencode.v34lift37.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz
wget -O data/gtf/gencode.v34.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz
gunzip data/annovar/humandb/hg38_phastConsElements100way.txt.gz
gunzip data/annovar/humandb/hg38_refGene.txt.gz
gunzip data/annovar/humandb/hg38_refGeneMrna.fa.gz
gunzip data/annovar/humandb/hg38_refGeneVersion.txt.gz
gunzip data/annovar/humandb/hg38_dbnsfp33a.txt.gz
gunzip data/annovar/humandb/hg38_dbnsfp33a.txt.idx.gz
gunzip data/annovar/humandb/hg38_dbnsfp41a.txt.gz
gunzip data/annovar/humandb/hg38_dbnsfp41a.txt.idx.gz
gunzip data/annovar/humandb/hg38_dbnsfp42a.txt.gz
gunzip data/annovar/humandb/hg38_gnomad30_genome.txt.gz
gunzip data/annovar/humandb/hg38_gnomad30_genome.txt.idx.gz
gunzip data/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip data/Homo_sapiens/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip data/gtf/gencode.v47lift37.annotation.gtf.gz
gunzip data/gtf/gencode.v47.annotation.gtf.gz
#Çë°²×°samtools£ºapt install samtools
samtools faidx data/Homo_sapiens/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa
samtools faidx data/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
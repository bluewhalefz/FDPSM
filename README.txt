## FDPSM framework
### Basic requirements
1. Python v3.9.12  
   LightGBM v4.4.0
   OpenFE v0.0.12

2. Use `conda env create -f fdpsm.yml` to create the conda environment is recommanded. 

3. SpliceAI(Use `conda env create -f spliceai.yml` to create the conda environment is recommanded).

4. AnnoVar (register required).


### Usage

1. Features derived from CADD (https://cadd.gs.washington.edu/) (GRCh38-v1.7)
   Note: After extracting the features to CADD, they need to be processed. Because the annotations obtained through CADD, some mutations may correspond to multiple lines of annotations.
   python process_CADDhg38.py --inputfile inputfile --cadd_file cadd_file --outfile outfile
   inputfile: Original file including basic information about the mutation.
   cadd_file: Files in csv format obtained through CADD annotations.
   outfile: Output file in csv format.

2. Features drived from SpliceAI, ANNOVAR and ChromHMM
   conda activate fdpsm
   python code/annotation_3tools.py --inputfile inputfile --outfile outfile
   
   Input format of input file
   |  Chr  | Start |  End  |  Ref  |  Alt  | 
   | ----- | ----- | ----- | ----- | ----- | 
   ÈçÎÄ¼þtest.csv
   
3. Features drived from MMsplice
   conda activate fdpsm
   ./features_MMSplice.sh  $1  $2  $3  $4  $5
   $1: Input file in vcf format.
   $2: hg38/hg19.
   $3: Output folder for output files.
   $4: The original file (in csv format) of your input file, including labels and other basic information.
   $5: Final output file with basic information about the mutation and the features predicted by MMSplice.
   
4. Make prediction
   python code/FDPSM_predict.py --train_file data/train8502_hg38_111fea_processed.csv --test_file data/test816_hg38_111fea_processed.csv --openfe_features model/openFE_train8502_111fea.features --selectedFea model/train8502_111fea_selected.csv --model_file model/train8502_111fea_FDPSM.model --outfile model/test816_111fea_FDPSM.csv
   
5. Model training
python code/FDPSM.py --train_file  data/train8502_hg38_111fea_processed.csv --selectedFea model_new/train8502_111fea_selected.csv --openfe_features model_new/openFE_train8502_111fea.features --model_file model_new/train8502_111fea_FDPSM.model --importance_file model_new/train8502_111fea_importance.csv --feaName_file model_new/train8501_111fea_feaName.txt  

   
   
   

   
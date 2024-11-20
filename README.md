## FDPSM framework
### Basic requirements
1. Use `conda env create -f fdpsm.yml` to create the conda environment is recommanded.  <br />
   (python &nbsp; v3.9.12 &emsp; LightGBM &nbsp; v4.4.0  &emsp;   OpenFE &nbsp; v0.0.12)<br />

2. SpliceAI(Use `conda env create -f spliceai.yml` to create the conda environment is recommanded).<br />

3. mmsplice<br />
   conda create -n mmsplice python=3.8<br />
   pip install cyvcf2 cython<br />
   pip install mmsplice<br />
   pip install numpy==1.26.4 <br />

4. Download file: bash download.sh<br />

### Usage

1. Features derived from CADD (https://cadd.gs.washington.edu/) (GRCh38-v1.7)<br />
   Note: After extracting the features to CADD, they need to be processed. Because the annotations obtained through CADD, some mutations may correspond to multiple lines of annotations.<br />
   python process_CADDhg38.py --inputfile inputfile --source_file source_file --outfile outfile<br />
   parameter description:<br />
      inputfile: Files in csv format obtained through CADD annotations. <br />(Feature files obtained from cadd website annotations, manual conversion of tsv format to csv format)<br />
      source_file: Original file including basic information about the mutation<br />
      outfile: Output file in csv format.<br />
   example:  python code/process_CADDhg38.py --inputfile data/example/test_CADD_hg38.csv --source_file data/example/test.csv --outfile example/test_CADD.csv<br />
   
2. Features drived from SpliceAI, ANNOVAR and ChromHMM<br />
   
   Annotate mutations: bash annotation_3tools.sh  $1 $2 $3<br />
   $1: the path of the input file (data/example/test.csv), <br />
   $2: the path of the output file (data/example/test_out.csv) <br />
   $3: the name of the input file (test)(no contain extensions).<br />
   
   Input format of input file<br />
   |  Chr  | Start |  End  |  Ref  | Alt | 
   |:----|:----|:----|:----|:----|
   | chr15 | 89260843 | 89260843 | G | A |
   |...|...|...|...|...| <br />
   "Chr" indicates the chromosome number, "Start" indicates the specific location where the mutation begins, <br />
   "End" indicates the specific location where the mutation ends, "Ref" indicates the base before the mutation, <br />
   and "Alt" indicates the base after the mutation<br />
   
   example: bash annotation_3tools.sh example/test.csv example/test_out.csv test<br />
   
3. Features drived from MMsplice<br />
   conda activate mmsplice<br />
   ./features_MMSplice.sh  $1  $2  $3  $4  $5<br />
   parameter description:<br />
      $1: Input file in vcf format.<br />
      $2: hg38/hg19.<br />
      $3: Output folder for output files.<br />
      $4: The original file (in csv format) of your input file, including labels and other basic information.<br />
      $5: Final output file with basic information about the mutation and the features predicted by MMSplice.<br />
   example: bash features_MMSplice.sh data/example/test_hg38.vcf hg38 data/example/out data/example/out/test_out.csv data/example/out/test_mmsplice.csv<br />
   Note: "$4" can directly use the output file from Usage2 <br />
        View log files when you encounter errors during feature generation -> test.std.avinput.features-mmsplice.tsv.log<br />
        (Note: Whether the file has the executable permission)<br />

4. Merge features<br />
   Combine all features into a single file (csv) that includes the mutation labels and other information, placing the non-feature information columns in front of the feature columns.<br />
   Also, sort the feature columns by the fea_list list (this step is necessary).<br />
   fea_list = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'MMS_delta_logit_psi', 'MMS_ref_acceptorIntron', 'MMS_alt_acceptorIntron', 'MMS_ref_acceptor', 'MMS_alt_acceptor', 'MMS_ref_exon', 'MMS_alt_exon', 'MMS_ref_donor', 'MMS_alt_donor', 'MMS_ref_donorIntron', 'MMS_alt_donorIntron', 'MMS_pathogenicity', 'MMS_efficiency', 'Dist2Mutation', 'Freq100bp', 'Rare100bp', 'Sngl100bp', 'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp', 'Rare10000bp', 'Sngl10000bp', 'priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP', 'verPhyloP', 'bStatistic', 'ZooPriPhyloP', 'ZooVerPhyloP', 'GC', 'CpG', 'minDistTSS', 'minDistTSE', 'GerpN', 'GerpS', 'RegSeq0', 'RegSeq1', 'RegSeq2', 'RegSeq3', 'RegSeq4', 'RegSeq5', 'RegSeq6', 'RegSeq7', 'gdi', 'gdi_phred', 'cHmm_E1', 'cHmm_E2', 'cHmm_E3', 'cHmm_E4', 'cHmm_E5', 'cHmm_E6', 'cHmm_E7', 'cHmm_E8', 'cHmm_E9', 'cHmm_E10', 'cHmm_E11', 'cHmm_E12', 'cHmm_E13', 'cHmm_E14', 'cHmm_E15', 'cHmm_E16', 'cHmm_E17', 'cHmm_E18', 'cHmm_E19', 'cHmm_E20', 'cHmm_E21', 'cHmm_E22', 'cHmm_E23', 'cHmm_E24', 'cHmm_E25', 'EncodeH3K4me1_sum', 'EncodeH3K4me1_max', 'EncodeH3K4me2_sum', 'EncodeH3K4me2_max', 'EncodeH3K4me3_sum', 'EncodeH3K4me3_max', 'EncodeH3K9ac_sum', 'EncodeH3K9ac_max', 'EncodeH3K9me3_sum', 'EncodeH3K9me3_max', 'EncodeH3K27ac_sum', 'EncodeH3K27ac_max', 'EncodeH3K27me3_sum', 'EncodeH3K27me3_max', 'EncodeH3K36me3_sum', 'EncodeH3K36me3_max', 'EncodeH3K79me2_sum', 'EncodeH3K79me2_max', 'EncodeH4K20me1_sum', 'EncodeH4K20me1_max', 'EncodeH2AFZ_sum', 'EncodeH2AFZ_max', 'EncodeDNase_sum', 'EncodeDNase_max', 'Hmec', 'Hsmm', 'Huvec', 'K562', 'Nhek', 'Nhlf']<br />
   Refer to the code below:<br />
   '''<br />
      colnames = df.columns.tolist()
      colnames = [name.replace("-", "_") for name in colnames]
      df_match.columns = colnames
      df_fea = df[fea_list].reindex(columns = fea_list)<br />
      df_match = pd.concat([df.iloc[:, :10], df_fea], axis = 1)  # "df.iloc[:, :10]" indicates the non-feature information column.<br />
      df_match.to_csv(outfile, index=False)<br />
   '''

5. Make prediction<br />
   python code/FDPSM_predict.py --train_file data/train8502_hg38_111fea_processed.csv --test_file data/test816_hg38_111fea_processed.csv --openfe_features model/openFE_train8502_111fea.features --selectedFea model/train8502_111fea_selected.csv --model_file model/train8502_111fea_FDPSM.model --outfile model/test816_111fea_FDPSM.csv<br />
   Note: The non-feature information column of the training set and the test set has 10 columns. If it is not 10 columns, please pay attention to modify the corresponding part of the code.<br />

6. Model training<br />
   python code/FDPSM.py --train_file  data/train8502_hg38_111fea_processed.csv --selectedFea model_new/train8502_111fea_selected.csv --openfe_features model_new/openFE_train8502_111fea.features --model_file model_new/train8502_111fea_FDPSM.model --importance_file model_new/train8502_111fea_importance.csv --feaName_file model_new/train8501_111fea_feaName.txt  <br />
   Note: The non-feature information column of the training set and the test set has 10 columns. If it is not 10 columns, please pay attention to modify the corresponding part of the code. <br />
   
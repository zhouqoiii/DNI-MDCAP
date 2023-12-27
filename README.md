# DNI-MDCAP
DNI-MDCAP: Improvement of Causal MiRNA-Disease Association Prediction Based on Deep Network Imputation.
## Dependencies
+ python = 3.9
+ sklearn
+ networkx
+ R = 4.2.1
## Instructions
+ '4models_comparison' folder contains all prediction scores of different models based on filtered training and testing sets of DNI-MDCAP, in which label 1, 2, 3 represent causal, non-causal, non-disease respectively.
+ 'data' folder contains all datasets for DNI-MDCAP.
+ 'ge' folder is cited from: Wu XB, Zhou Y: GE-Impute: graph embedding-based imputation for single-cell RNA-seq data. Brief Bioinform 2022, 23(5):bbac313.
+ 'network' folder contains the miRNA metrics and miRNA-disease metrics after network imputation.
+ 'result' folder contains all prediction results of DNI-MDCAP.
## Tutorial
```
#Step1 load the similarity matrix of miRNA, disease, and known causal miRNA-disease association in 'data' folder.
rawfile=pd.read_csv("input_file.txt",sep="\t",index_col=None,header=None) #miR_sim_matrix, disease_sim_matrix, miR_disease_matrix

#Step2 Network Imputation in Inputation Network.ipynb
embeddings_df=graph_embedding(rawfile) 
result=imputation(embeddings_df,rawfile) 

#Step3 Semi-supervised learning framework in ie_l1norm.R
miR_info <- LP_norm1(miR_sim_matrix, miR_disease_matrix, con_condition)
disease_info <- LP_norm1(disease_sim_matrix, t(miR_disease_matrix), con_condition)
prediction <- miR_info$Q * factor + t(disease_info$Q) * (1 - factor)

```

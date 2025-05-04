# scIHDB
ScIHDB is a comprehensive resource for ischemic heart disease (IHD) research with an well-annotated single-cell atlas and a genetic, pharmacological, and comorbid knowledge graph.
![cover letter](https://github.com/user-attachments/assets/4f64639a-e68c-494d-b21f-6da1f05f9421)


# Introduction to scIHDB
**S**ingle-**C**ell **I**ntegrated **D**ata**B**ase of **I**schemic **H**eart Disease (**scIHDB**, https://xomics.com.cn/ihdb) is a comprehensive resource for ischemic heart disease that allows spatial and temporal interpretations. ScIHDB is composed of a single-cell atlas of IHD (**IHDAtlas**) encompassing over 1.4 million cells that can be categorized into 40 common cell types, and a IHD knowledge graph (**IHDKG**) with manually curated 6,035 CHD-related genes, 3,157 drugs, 80 comorbidities, and gene-disease associations and molecular interactions between them. To make scIHDB a reliable reference for IHD research, we performed collective meta-program (MP) analysis and generated 16 MPs associated with different murine disease models. Each MP is a representitive of a group of gene programs under disease states. Besides, a pipeline was constructed from a series of state-of-the-art computational tools to analyze each dataset individually. All analysis results, including differential gene expression, cell-cell interaction, pathway activeness, were visulized and publicly available on our website. Our IHDKG was so far the largest manually curated knowledge graph for ischemic heart disease. It provides the most genes, drugs, and gene-disease associations with supporting literature comparing to other dedicated databases. Besides, IHDKG is a "dynamic" knowledge graph with each gene been annotated with its significance in cell types under different disease conditions. Our website provides a network creation and manipulation tool for users to build subgraphs for download and further analysis with tools such as networkx and SNAP.

To run our analysis pipeline, we suggest to install Anaconda on Windows PC and create a new environment using the sc.yml configuration to install all related packages.
```
conda env create -f sc.yml -n sc
conda activate sc
```
Besides, an R environment should be configured to include all related packages. The required R packages include:  
```
NMF==0.27  
scry==1.16.0  
scran==1.32.0  
Biobase==2.62.0  
BiocParallel==1.36.0  
DESeq2==1.38.3  
anndata==0.7.5.6  
Matrix==1.6-5  
zinbwave==1.28.0  
caret==7.0-1  
msigdbr==7.5.1  
gsdensity==0.1.3  
```
Before preprocessing scRNA-seq data, users should reformat and organize their scRNA-seq matrix data into one of the following structures:  
```
1. 10x  
   GSEXXXXXX  
    ├─GSMXXXXXXX  
    │      barcodes.tsv.gz  
    │      features.tsv.gz  
    │      matrix.mtx.gz  
    │  
    ├─GSMXXXXXXX  
    │      barcodes.tsv.gz  
    │      features.tsv.gz  
    │      matrix.mtx.gz  
    ...  
2. 10x tiled  
   GSEXXXXXX  
       GSMXXXXXXX_barcodes.tsv.gz  
       GSMXXXXXXX_features.tsv.gz  
       GSMXXXXXXX_matrix.mtx.gz  
       GSMXXXXXXX_barcodes.tsv.gz  
       GSMXXXXXXX_features.tsv.gz  
       GSMXXXXXXX_matrix.mtx.gz  
       ...  
3. csv  
   GSEXXXXXX  
        GSMXXXXXXX.csv  
        GSMXXXXXXX.csv  
        ...  
5. txt  
   GSEXXXXXX  
        GSMXXXXXXX.txt  
        GSMXXXXXXX.txt  
        ...  
6. h5  
   GSEXXXXXX  
        GSMXXXXXXX.h5  
        GSMXXXXXXX.h5  
        ...  
7. h5ad  
   GSEXXXXXX  
        GSMXXXXXXX.h5ad  
        GSMXXXXXXX.h5ad  
        ...  
```
Put all dataset folders under the same folder. Run "preprocess_QC.ipynb".  

For single-cell analysis of individual samples, run "Annotation.py" first for every preprocessed files, then run "HVG.R" for all files, finally "PATHWAY.R", "NMF.R", and "run_EMBEDandCCI.py".  
For meta-program analysis, run "generate_robust_programs.ipynb" first, then run "generate_meta_programs.R", which is modified from the orginal "Generate_Meta_Programs.R".  
As we only provide the modified code for meta-program generation, for complete and orginal codes please refer to https://github.com/tiroshlab/3ca/tree/main/ITH_hallmarks

# About
Should you have any questions, please feel free to contact the auther, Mr. Tianhao Wang (woloorn@zju.edu.cn)

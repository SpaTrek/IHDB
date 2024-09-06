# idCHD
**I**ntegrated **D**atabase of **C**oronary **H**eart **D**isease (**idCHD**, https://xomics.com.cn/idchd) is a comprehensive resource and online toolkit for coronary heart disease-centric medical research. It curates CHD-related studies and single-cell datasets published over decades and provides interactive visualization of analytical results, with a tool for multi-level interaction network construction. The database encompasses 585,412 cells categorized into 17 major classes or 126 subclasses from high-quality scRNA-seq and snRNA-seq datasets. Differential gene expression, cell-cell interaction, gene set enrichment analysis, gene program factorization, and meta-analysis of gene programs revealed a single-cell landscape of disease progression. Additionally, with manually curated 6,035 CHD-related genes, 3,157 drugs, and ~150,000 molecular interactions and gene-comorbidity associations, multi-level networks can be constructed and downloaded for further analysis. In summary, idCHD serves as a multi-dimensional resource for CHD-centric medical research, provides deep insights into genetic and cellular architecture of complex diseases.

![cover letter fig 1](https://github.com/user-attachments/assets/73e19c35-5de9-4bb6-9bf5-4cf90350b000)

To run the single-cell RNA-seq pipeline, we suggest to install Anaconda on a Windows PC and create a new environment using the sc.yml configuration to install all related packages.
```
conda env create -f sc.yml -n sc
conda activate sc
```
Besides, an R environment should be configured to include all related packages. The required R packages include:
NMF==0.27
scry==1.16.0
scran==1.32.0
Biobase==2.62.0
iocParallel==1.36.0

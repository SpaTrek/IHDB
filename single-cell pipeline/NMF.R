library(scry); library(Biobase); library(NMF); library(BiocParallel)

library(anndata)
library(Matrix)

library(reticulate)
use_python("/path/to/", required = T)
py_config()

path <- "/path/to/all/preprocessed/h5ad/files/"
files <- dir(path)

human_datasets<-c()#list of human dataset IDs


for (file in files) {
  print(paste("Start factorizng ",file,sep = ""))
  Sample <- substr(file,1,10)
    if (TRUE) {
      infile = paste(path,paste(file,".h5ad",sep=""),sep = "")
      g <- read_h5ad(infile)
      g_X_t <- Matrix::t(g$layers["lognormal"])
      g_X_t <- g_X_t[!(startsWith(rownames(g_X_t), "ERCC-") | (rownames(g_X_t) %in% c("EGFP","YFP","DsRed","tdTomato","ZsGreen1"))),]
      
      sce <- devianceFeatureSelection(g_X_t)
      highly_variable_matrix <- g_X_t[order(sce, decreasing = T)[1:5000],]
      highly_variable_matrix <- apply(highly_variable_matrix, 1, function(x)x-(mean(x)))
      highly_variable_matrix[highly_variable_matrix<0] <- 0
      res <- nmf(highly_variable_matrix, 30, method="snmf/r", seed="nndsvd")
      weight <- t(coef(res))
      embeddings <- basis(res)
      colnames(weight) <- paste(rep(paste(Sample,"_rank4_9_nruns10.RDS.30.",sep=""),30), as.character(1:30),sep="")
      colnames(embeddings) <- paste(rep("nmf_",30), as.character(1:30),sep="")
      write.csv(weight, paste("/path/to/NMF/folder/",Sample,"_weights.csv", sep=""))
      write.csv(embeddings, paste("/path/to/NMF/folder/",Sample,"_module.csv", sep=""))
    }
#  }
}

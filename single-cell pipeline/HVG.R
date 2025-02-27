library(scry)
library(Seurat)
library(DESeq2)
library(edgeR)
library(anndata)
library(Matrix)
library(zinbwave)
library(reticulate)
library(dplyr)
library(caret)
use_python("/path/to/python/binary", required = T)
py_config()

path <- "/path/to/all/preprocessed/h5ad/files/"
files <- dir(path)

human_datasets<-c()#list of human dataset IDs

for (file in files) {
  if (!is.na(file)) {
    print(paste("Start analyzing ",file,sep = ""))
    Sample <- substr(file,1,10)
    if (!(paste(Sample,"_hv.csv",sep="") %in% dir(path))) {
      infile <- paste(path,file,sep = "")
      g <- read_h5ad(infile)
      g_X_t <- Matrix::t(g$layers["counts"])
      g_X_t <- round(g_X_t,0)
      g_X_t <- g_X_t[!(startsWith(rownames(g_X_t), "ERCC-") | (rownames(g_X_t) %in% c("EGFP","YFP","DsRed","tdTomato","ZsGreen1"))),]
      g_X_t <- g_X_t[(rowSums(g_X_t)>50), ]
      celltype <- g$obs$celltypist_cell_label
      if ((length(levels(celltype))==1) & (dim(g_X_t)[2]<=150)) {
        print("Small sample. Skip.")
        next
      } 
      else {
        print("Sample contains >150 cells. OK.")
        sce <- devianceFeatureSelection(g_X_t)
        
        seu <- CreateSeuratObject(g_X_t)
        seu$celltype <- celltype
        seu <- seu %>% NormalizeData %>% FindVariableFeatures(nfeatures = 5000, selection.method = "vst") %>% ScaleData
        seu <- RunPCA(seu, features = VariableFeatures(object = seu)) %>% FindNeighbors(dims = 1:30) %>% RunUMAP(dims = 1:30)
        seu@reductions$umap@cell.embeddings <- g$obsm$X_umap
        rownames(seu@reductions$umap@cell.embeddings) <- colnames(seu)
        colnames(seu@reductions$umap@cell.embeddings) <- c("umap_1","umap_2")
        seu$subtype <- as.character(celltype)
        
        sava <- data.frame(baseMean=0, log2FoldChange=0, lfcSE=0, pvalue=0, padj=0, cell_type=0)
        for (ct in levels(celltype)) {
          print(paste("Analyzing ", ct, ".", sep = ""))
          
          if (table(seu$celltype)[ct] > 150) {
            print("Celltype with >150 cells detected. Conduct subtyping...")
            sbst <- seu[ ,seu$celltype==ct]
            
            sbst <- FindClusters(sbst, resolution = 0.3)
            if (length(levels(sbst$seurat_clusters))==1) {
              print("Homogenous celltype. Conduct celltype-level regression.")
              seu@meta.data[colnames(sbst), ]$subtype <- as.character(ct)
              
              df <- data.frame(celltype = as.character(seu$celltype))
              df$celltype[df$celltype!=as.character(ct)] <- 0
              df$celltype[df$celltype==as.character(ct)] <- 1
              
              obj <- SummarizedExperiment(assays=list(counts=g_X_t), colData=df)
              vars <- assay(obj) %>% log1p %>% rowVars
              names(vars) <- rownames(obj)
              vars <- sort(vars, decreasing = TRUE)
              vars <- names(vars)[1:5000]
              vars <- vars[!is.na(vars)]
              obj <- obj[vars,]
              assayNames(obj)[1] <- "counts"
              
              print("zinbFit...")
              zinb <- zinbFit(obj, K=2, epsilon=1e12)
              obj_zinb <- zinbwave(obj, fitted_model = zinb, K = 2, epsilon=1e12, observationalWeights = TRUE)
              counts(obj_zinb) <- as.matrix(counts(obj_zinb))
              counts(obj_zinb) <- counts(obj_zinb)+matrix(1, nrow = dim(counts(obj_zinb))[1], ncol = dim(counts(obj_zinb))[2])
              
              test <- try({
                print("DESeq2 regression...")
                dds <- DESeqDataSet(obj_zinb, design = ~ celltype)
                dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
                res <- lfcShrink(dds, contrast=c("celltype", 1, 0), type = "normal")
                result <- as.data.frame(res@listData)
                rownames(result) <- rownames(res)
                result <- result %>% arrange(desc(log2FoldChange), padj)
                savv <- head(result, 10)
                if ("stat" %in% colnames(savv)) {
                  savv <- savv[ ,-4]
                }
                savv$cell_type <- as.character(ct)
                
                sav_if <- intersect(rownames(savv), rownames(sava))
                if (length(sav_if)==0) {
                  sava <- rbind(sava, savv)
                } else {
                  for (rep_gene in sav_if) {
                    sava[rep_gene, ]$cell_type <- paste(sava[rep_gene, ]$cell_type, "|", savv[rep_gene, ]$cell_type)
                  }
                  savv <- savv[!(rownames(savv) %in% sav_if), ]
                  sava <- rbind(sava, savv)
                }
              })
              if (!('try-error' %in% class(test))) {
                print(paste("Error. Skip ",ct,".",sep = ""))
                next
              }
            }
            else {
              print("Multiple subtypes detected. Conduct subtype-level regression.")
              
              seu@meta.data[colnames(sbst), ]$subtype <- paste(rep(as.character(ct), dim(sbst)[2]), as.character(sbst$seurat_clusters), sep =".")
              
              df <- data.frame(celltype = sbst$seurat_clusters)
              dummy <- caret::dummyVars("~.", data = df)
              final_df <- data.frame(predict(dummy, newdata=df))
              final_df <- as.data.frame(lapply(final_df, function(x) factor(x)), row.names = rownames(final_df))
              
              obj <- SummarizedExperiment(assays=list(counts=g_X_t[ ,colnames(sbst)]),
                                          colData=final_df)
              subtypes <- levels(sbst$seurat_clusters)
              vars <- assay(obj) %>% log1p %>% rowVars
              names(vars) <- rownames(obj)
              vars <- sort(vars, decreasing = TRUE)
              vars <- names(vars)[1:3000]
              vars <- vars[!is.na(vars)]
              obj <- obj[vars,]
              assayNames(obj)[1] <- "counts"
              
              print("zinbFit...")
              zinb <- zinbFit(obj, K=2, epsilon=1e12)
              obj_zinb <- zinbwave(obj, fitted_model = zinb, K = 2, epsilon=1e12, observationalWeights = TRUE)
              counts(obj_zinb) <- as.matrix(counts(obj_zinb))
              counts(obj_zinb) <- counts(obj_zinb)+matrix(1, nrow = dim(counts(obj_zinb))[1], ncol = dim(counts(obj_zinb))[2])
              
              sava_backup <- sava
              
              test <- try({
                for (subtype in subtypes) {
                  print(paste("Regressing ", subtype, ".", sep = ""))
                  
                  dds <- DESeqDataSet(obj_zinb, design = as.formula(paste("~ celltype", subtype, sep = ".")))
                  
                  print("DESeq2 regression...")
                  dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
                  
                  res <- lfcShrink(dds, contrast=c(paste("celltype.",subtype,sep = ""), 1, 0), type = "normal")
                  
                  result <- as.data.frame(res@listData)
                  rownames(result) <- rownames(res)
                  result <- result %>% arrange(desc(log2FoldChange), padj)
                  savv <- head(result, 10)
                  if ("stat" %in% colnames(savv)) {
                    savv <- savv[ ,-4]
                  }
                  savv$cell_type <- paste(ct, subtype, sep = ".")
                  
                  sav_if <- intersect(rownames(savv), rownames(sava))
                  if (length(sav_if)==0) {
                    sava <- rbind(sava, savv)
                  } else {
                    for (rep_gene in sav_if) {
                      sava[rep_gene, ]$cell_type <- paste(sava[rep_gene, ]$cell_type, "|", savv[rep_gene, ]$cell_type)
                    }
                    savv <- savv[!(rownames(savv) %in% sav_if), ]
                    sava <- rbind(sava, savv)
                  }
                }
              })
              if (!('try-error' %in% class(test))) {
                print(paste("Done: ",ct,sep = ""))
                next
              } # Else, clustering resolution = 0.2
              else {
                print(paste("Error: subtype.",subtype," Switch to reso=0.2", sep = ""))
                
                sava <- sava_backup
                
                sbst <- FindClusters(sbst, resolution = 0.2)
                if (length(levels(sbst$seurat_clusters))==1) {
                  print("Homogenous celltype. Conduct celltype-level regression.")
                  seu@meta.data[colnames(sbst), ]$subtype <- as.character(ct)
                  
                  df <- data.frame(celltype = as.character(seu$celltype))
                  df$celltype[df$celltype!=as.character(ct)] <- 0
                  df$celltype[df$celltype==as.character(ct)] <- 1
                  
                  obj <- SummarizedExperiment(assays=list(counts=g_X_t), colData=df)
                  vars <- assay(obj) %>% log1p %>% rowVars
                  names(vars) <- rownames(obj)
                  vars <- sort(vars, decreasing = TRUE)
                  vars <- names(vars)[1:5000]
                  vars <- vars[!is.na(vars)]
                  obj <- obj[vars,]
                  assayNames(obj)[1] <- "counts"
                  
                  print("zinbFit...")
                  zinb <- zinbFit(obj, K=2, epsilon=1e12)
                  obj_zinb <- zinbwave(obj, fitted_model = zinb, K = 2, epsilon=1e12, observationalWeights = TRUE)
                  counts(obj_zinb) <- as.matrix(counts(obj_zinb))
                  counts(obj_zinb) <- counts(obj_zinb)+matrix(1, nrow = dim(counts(obj_zinb))[1], ncol = dim(counts(obj_zinb))[2])
                  
                  test <- try({
                    print("DESeq2 regression...")
                    dds <- DESeqDataSet(obj_zinb, design = ~ celltype)
                    dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
                    res <- lfcShrink(dds, contrast=c("celltype", 1, 0), type = "normal")
                    result <- as.data.frame(res@listData)
                    rownames(result) <- rownames(res)
                    result <- result %>% arrange(desc(log2FoldChange), padj)
                    savv <- head(result, 10)
                    if ("stat" %in% colnames(savv)) {
                      savv <- savv[ ,-4]
                    }
                    savv$cell_type <- as.character(ct)
                    
                    sav_if <- intersect(rownames(savv), rownames(sava))
                    if (length(sav_if)==0) {
                      sava <- rbind(sava, savv)
                    } else {
                      for (rep_gene in sav_if) {
                        sava[rep_gene, ]$cell_type <- paste(sava[rep_gene, ]$cell_type, "|", savv[rep_gene, ]$cell_type)
                      }
                      savv <- savv[!(rownames(savv) %in% sav_if), ]
                      sava <- rbind(sava, savv)
                    }
                  })
                  if (!('try-error' %in% class(test))) {
                    print(paste("Error. Skip ",ct,".",sep = ""))
                    next
                  }
                }
                else {
                  print("Multiple subtypes detected. Conduct subtype-level regression.")
                  
                  seu@meta.data[colnames(sbst), ]$subtype <- paste(rep(as.character(ct), dim(sbst)[2]), as.character(sbst$seurat_clusters), sep =".")
                  
                  df <- data.frame(celltype = sbst$seurat_clusters)
                  dummy <- caret::dummyVars("~.", data = df)
                  final_df <- data.frame(predict(dummy, newdata=df))
                  final_df <- as.data.frame(lapply(final_df, function(x) factor(x)), row.names = rownames(final_df))
                  
                  obj <- SummarizedExperiment(assays=list(counts=g_X_t[ ,colnames(sbst)]),
                                              colData=final_df)
                  subtypes <- levels(sbst$seurat_clusters)
                  vars <- assay(obj) %>% log1p %>% rowVars
                  names(vars) <- rownames(obj)
                  vars <- sort(vars, decreasing = TRUE)
                  vars <- names(vars)[1:3000]
                  vars <- vars[!is.na(vars)]
                  obj <- obj[vars,]
                  assayNames(obj)[1] <- "counts"
                  
                  print("zinbFit...")
                  zinb <- zinbFit(obj, K=2, epsilon=1e12)
                  obj_zinb <- zinbwave(obj, fitted_model = zinb, K = 2, epsilon=1e12, observationalWeights = TRUE)
                  counts(obj_zinb) <- as.matrix(counts(obj_zinb))
                  counts(obj_zinb) <- counts(obj_zinb)+matrix(1, nrow = dim(counts(obj_zinb))[1], ncol = dim(counts(obj_zinb))[2])
                  
                  sava_backup <- sava
                  
                  test <- try({
                    for (subtype in subtypes) {
                      print(paste("Regressing ", subtype, ".", sep = ""))
                      
                      dds <- DESeqDataSet(obj_zinb, design = as.formula(paste("~ celltype", subtype, sep = ".")))
                      
                      print("DESeq2 regression...")
                      dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
                      
                      res <- lfcShrink(dds, contrast=c(paste("celltype.",subtype,sep = ""), 1, 0), type = "normal")
                      
                      result <- as.data.frame(res@listData)
                      rownames(result) <- rownames(res)
                      result <- result %>% arrange(desc(log2FoldChange), padj)
                      savv <- head(result, 10)
                      if ("stat" %in% colnames(savv)) {
                        savv <- savv[ ,-4]
                      }
                      savv$cell_type <- paste(ct, subtype, sep = ".")
                      
                      sav_if <- intersect(rownames(savv), rownames(sava))
                      if (length(sav_if)==0) {
                        sava <- rbind(sava, savv)
                      } else {
                        for (rep_gene in sav_if) {
                          sava[rep_gene, ]$cell_type <- paste(sava[rep_gene, ]$cell_type, "|", savv[rep_gene, ]$cell_type)
                        }
                        savv <- savv[!(rownames(savv) %in% sav_if), ]
                        sava <- rbind(sava, savv)
                      }
                    }
                  })
                  if (!('try-error' %in% class(test))) {
                    print(paste("Done: ",ct,sep = ""))
                    next
                  } # Else, do celltype-level
                  else {
                    print(paste("Error: ",subtype,"! ","Switch to homo", sep = ""))
                    if (length(levels(celltype))==1) {
                      print("Homo sample with >150 cells. Skip")
                    }
                    else {
                      sava <- sava_backup
                      
                      print("Homogenous celltype. Conduct celltype-level regression.")
                      seu@meta.data[colnames(sbst), ]$subtype <- as.character(ct)
                      
                      df <- data.frame(celltype = as.character(seu$celltype))
                      df$celltype[df$celltype!=as.character(ct)] <- 0
                      df$celltype[df$celltype==as.character(ct)] <- 1
                      
                      obj <- SummarizedExperiment(assays=list(counts=g_X_t), colData=df)
                      vars <- assay(obj) %>% log1p %>% rowVars
                      names(vars) <- rownames(obj)
                      vars <- sort(vars, decreasing = TRUE)
                      vars <- names(vars)[1:5000]
                      vars <- vars[!is.na(vars)]
                      obj <- obj[vars,]
                      assayNames(obj)[1] <- "counts"
                      
                      print("zinbFit...")
                      zinb <- zinbFit(obj, K=2, epsilon=1e12)
                      obj_zinb <- zinbwave(obj, fitted_model = zinb, K = 2, epsilon=1e12, observationalWeights = TRUE)
                      counts(obj_zinb) <- as.matrix(counts(obj_zinb))
                      counts(obj_zinb) <- counts(obj_zinb)+matrix(1, nrow = dim(counts(obj_zinb))[1], ncol = dim(counts(obj_zinb))[2])
                      
                      test <- try({
                        print("DESeq2 regression...")
                        dds <- DESeqDataSet(obj_zinb, design = ~ celltype)
                        dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
                        res <- lfcShrink(dds, contrast=c("celltype", 1, 0), type = "normal")
                        result <- as.data.frame(res@listData)
                        rownames(result) <- rownames(res)
                        result <- result %>% arrange(desc(log2FoldChange), padj)
                        savv <- head(result, 10)
                        if ("stat" %in% colnames(savv)) {
                          savv <- savv[ ,-4]
                        }
                        savv$cell_type <- as.character(ct)
                        
                        sav_if <- intersect(rownames(savv), rownames(sava))
                        if (length(sav_if)==0) {
                          sava <- rbind(sava, savv)
                        } else {
                          for (rep_gene in sav_if) {
                            sava[rep_gene, ]$cell_type <- paste(sava[rep_gene, ]$cell_type, "|", savv[rep_gene, ]$cell_type)
                          }
                          savv <- savv[!(rownames(savv) %in% sav_if), ]
                          sava <- rbind(sava, savv)
                        }
                      })
                      if (!('try-error' %in% class(test))) {
                        print(paste("Error. Skip ",ct,".",sep = ""))
                        next
                      }
                    }
                  }
                }
              }
            }
          } else {
            print("Small celltype detected. Conduct celltype-level regression.")
            sbst <- seu[ ,seu$celltype==ct]
            seu@meta.data[colnames(sbst), ]$subtype <- as.character(ct)
            
            df <- data.frame(celltype = as.character(seu$celltype))
            df$celltype[df$celltype!=as.character(ct)] <- 0
            df$celltype[df$celltype==as.character(ct)] <- 1
            
            obj <- SummarizedExperiment(assays=list(counts=g_X_t), colData=df)
            vars <- assay(obj) %>% log1p %>% rowVars
            names(vars) <- rownames(obj)
            vars <- sort(vars, decreasing = TRUE)
            vars <- names(vars)[1:5000]
            vars <- vars[!is.na(vars)]
            obj <- obj[vars,]
            assayNames(obj)[1] <- "counts"
            
            print("zinbFit...")
            zinb <- zinbFit(obj, K=2, epsilon=1e12)
            obj_zinb <- zinbwave(obj, fitted_model = zinb, K = 2, epsilon=1e12, observationalWeights = TRUE)
            counts(obj_zinb) <- as.matrix(counts(obj_zinb))
            counts(obj_zinb) <- counts(obj_zinb)+matrix(1, nrow = dim(counts(obj_zinb))[1], ncol = dim(counts(obj_zinb))[2])
            
            test <- try({
              dds <- DESeqDataSet(obj_zinb, design = ~ celltype)
              
              print("DESeq2 regression...")
              dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
              res <- lfcShrink(dds, contrast=c("celltype", 1, 0), type = "normal")
              result <- as.data.frame(res@listData)
              rownames(result) <- rownames(res)
              result <- result %>% arrange(desc(log2FoldChange), padj)
              savv <- head(result, 10)
              if ("stat" %in% colnames(savv)) {
                savv <- savv[ ,-4]
              }
              savv$cell_type <- as.character(ct)
              
              sav_if <- intersect(rownames(savv), rownames(sava))
              if (length(sav_if)==0) {
                sava <- rbind(sava, savv)
              } else {
                for (rep_gene in sav_if) {
                  sava[rep_gene, ]$cell_type <- paste(sava[rep_gene, ]$cell_type, "|", savv[rep_gene, ]$cell_type)
                }
                savv <- savv[!(rownames(savv) %in% sav_if), ]
                sava <- rbind(sava, savv)
              }
            })
            if (!('try-error' %in% class(test))) {
              print(paste("Error. Skip ",ct,".",sep = ""))
              next
            }
          }
        }
        print("Analysis done! Saving files.")
        
        seu$subtype <- as.factor(seu$subtype)
        g$obs$subtype <- seu$subtype
        write_h5ad(g, paste(path, file, sep = ""))
        
        sava$gene <- rownames(sava)
        sava <- sava[-1, ]
        sava <- sava[ , c(7,6,2,5)]
        colnames(sava) <- c("gene","cell_type","logFC","pval_adj")
        sava$logFC <- round(sava$logFC, 2)
        
        write.csv(sava, paste(path,Sample,"_hv.csv", sep = ""), row.names = F, eol = "\r\n", quote = FALSE)
        
        
        mat <- t(seu@assays$RNA@layers$data)
        rownames(mat) <- colnames(seu)
        colnames(mat) <- rownames(seu)
        
        genehv <- mat[ ,intersect(sava$gene, colnames(mat))]
        genehv <- round(genehv, 2)
        genehv <- as.data.frame(genehv)
        genehv <- cbind("cell"=rownames(genehv), genehv)
        write.csv(genehv, paste(path,Sample,"_marker.csv", sep = ""), eol = "\r\n", row.names = F, quote = FALSE)
        
        geneapp <- mat[ ,c(intersect(sava$gene, colnames(mat)), intersect(names(sce)[order(sce, decreasing = T)[1:500]],colnames(mat)))]
        geneapp <- round(geneapp, 2)
        geneapp <- as.data.frame(geneapp)
        geneapp <- cbind("cell"=rownames(geneapp), geneapp)
        write.csv(geneapp, paste(path,Sample,"_app.csv", sep = ""), eol = "\r\n", row.names = F, quote = FALSE)
      }
    } else {
      print("Already analyzed!")
    }
  }
  gc()
}

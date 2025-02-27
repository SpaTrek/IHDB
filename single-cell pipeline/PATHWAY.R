library(Seurat)
library(scry)
library(anndata)
library(Matrix)
library(reticulate)
library(dplyr)
library(gsdensity)
library(msigdbr)
library(future)
library(future.apply)
use_python("/path/to/python/binary", required = T)
py_config()

path <- "/path/to/all/output/csv/files/"
files <- dir(path)
ot <- dir(path)

human_datasets<-c()#list of human dataset IDs

plan("multisession", workers=12)

for (file in files) {
  if (grepl("_id.csv",file)) {
    file <- substring(file,1,10)
    file <- paste(file,".h5ad",sep = "")
    if (!is.na(file) & length(which(ot==paste(substr(file,1,10),"_pathway.csv",sep = "")))==0) {
      Sample <- substr(file,1,10)
      if (TRUE) {
        print(paste("Start analyzing ",file,sep = ""))
      #if ((Sample %in% human_datasets)) {
        infile <- paste("/path/to/all/h5ad/files/",file,sep = "")
        g <- read_h5ad(infile)
        g_X_t <- Matrix::t(g$layers["counts"])
        g_X_t <- g_X_t[!(startsWith(rownames(g_X_t), "ERCC-") | (rownames(g_X_t) %in% c("EGFP","YFP","DsRed","tdTomato","ZsGreen1"))),]
        g_X_t <- g_X_t[(rowSums(g_X_t)>10), ]
        celltype <- g$obs$celltypist_cell_label
        sce <- devianceFeatureSelection(g_X_t)
        
        seu <- CreateSeuratObject(g_X_t)
        seu$subtype <- g$obs$subtype_revised
        ce <- compute.mca(object = seu)
        if ((Sample %in% human_datasets)) {
          gmt <- clusterProfiler::read.gmt("/pathway/to/human.gmt")
        } else {
          gmt <- clusterProfiler::read.gmt("/pathway/to/mouse.gmt")
        }
        gene.set.list <- list()
        for (gene.set.name in unique(gmt$term)) {
          gene.set.list[[gene.set.name]] <- gmt[gmt$term %in% gene.set.name, ]$gene
        }
        genes <- sapply(gene.set.list, function(x) paste(x, collapse = ", "))
        gene.set.list.df <- cbind(gene.set = names(gene.set.list), genes = genes)
        rownames(gene.set.list.df) <- 1:nrow(gene.set.list.df)
  
        res <- compute.kld(coembed = ce, 
                           genes.use = intersect(rownames(ce), rownames(seu)), 
                           n.grids = 100, 
                           gene.set.list = gene.set.list,
                           gene.set.cutoff = 3, 
                           n.times = 100
                           )
        
        gene.set.deviated <- res[res$p.adj < 0.001, ]$gene.set
        cells <- colnames(seu)
        el <- compute.nn.edges(coembed = ce, nn.use = 300)
        
        cv.df <- run.rwr.list(el = el, 
                              gene_set_list = gene.set.list[gene.set.deviated],
                              cells = cells)
        
        
        cv.df$subtype <- seu$subtype
        cv.all <- list()
        for (geneset in gene.set.deviated) {
          cv.all[[geneset]] <- aggregate(cv.df[[geneset]], by = list(type=cv.df$subtype), mean)$x
        }
        cv.all <- as.data.frame(cv.all)
        rownames(cv.all) <- aggregate(cv.df[[geneset]], by = list(type=cv.df$subtype), mean)$type
        cv.all <- cv.all[,!colSums(cv.all)==0]
        
        if (dim(cv.all)[1]!=1) {
          p <- pheatmap::pheatmap(as.matrix(cv.all), 
                             scale = "column", 
                             labels_col = F, 
                             clustering_distance_cols = "manhattan")
          
          cell_order <- p$tree_row$labels[p$tree_row$order]
          pathway_order <- p$tree_col$labels[p$tree_col$order]
          cv.all <- scale(cv.all)
          cv.all <- t(cv.all)
          cv.all <- cv.all[pathway_order,cell_order]
        } else {
          cv.all <- scale(cv.all)
          cv.all <- t(cv.all)
        }
        cv.all <- as.data.frame(cv.all)
        cv.all <- round(cv.all, 3)
        cv.all <- cbind("pathway"=rownames(cv.all), cv.all)
        
        write.csv(cv.all, paste(path,Sample,"_pathway.csv", sep = ""), row.names = F, quote = F, eol = "\r\n")
      }
    }
  }
  else {
    next
  }
}
      
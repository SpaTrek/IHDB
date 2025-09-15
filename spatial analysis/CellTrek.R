options(stringsAsFactors = F)
library(CellTrek)
library(dplyr)
library(anndata)
library(Matrix)
library(Seurat)
library(reticulate)
use_python("/slurm/home/yrd/liaolab/wangtianhao/anaconda3/envs/sc3.8.5/bin/python", required = T)
py_config()

g <- read_h5ad("~/data/idCHD/combined_cell2location.h5ad")
sc_ref <- CreateSeuratObject(t(g$X), meta.data = g$obs)
sc_ref <- RenameCells(sc_ref, new.names=make.names(Cells(sc_ref)))

files <- list.files("~/data/idCHD/spatial/")
for (file in files[47:54]) {
  g <- read_h5ad(paste("~/data/idCHD/spatial/",file,sep = ""))
  seu <- CreateSeuratObject(t(g$X), assay = "Spatial", meta.data = g$obs, )
  
  scale.factors <- list(fiducial_diameter_fullres = g$uns$spatial[[names(g$uns$spatial)]]$scalefactors$fiducial_diameter_fullres,
                        tissue_hires_scalef = g$uns$spatial[[names(g$uns$spatial)]]$scalefactors$tissue_hires_scalef,
                        spot_diameter_fullres = g$uns$spatial[[names(g$uns$spatial)]]$scalefactors$spot_diameter_fullres)
  
  
  
  tissue_positions_list <- data.frame(row.names = colnames(seu),
                                      tissue = 1,
                                      row = g$obsm$spatial[,2], col = g$obsm$spatial[,1],
                                      imagerow = g$obsm$spatial[,2], imagecol = g$obsm$spatial[,1])
  
  
  spatialObj <- new(Class = 'VisiumV1',
                    image = g$uns$spatial[[names(g$uns$spatial)]]$images$hires,
                    scale.factors = scalefactors(spot = scale.factors$spot_diameter_fullres,
                                                 fiducial = scale.factors$fiducial_diameter_fullres,
                                                 hires = scale.factors$tissue_hires_scalef,
                                                 lowres = scale.factors$tissue_hires_scalef),
                    coordinates = tissue_positions_list,
                    spot.radius = scale.factors$spot_diameter_fullres)
  
  
  spatialObj <- spatialObj[Cells(seu)]
  DefaultAssay(spatialObj) <- 'Spatial'
  seu[['slice1']] <- spatialObj
  
  
  ## Rename the cells/spots with syntactically valid names
  seu <- RenameCells(seu, new.names=make.names(Cells(seu)))
  
  traint <- CellTrek::traint(st_data=seu, sc_data=sc_ref, sc_assay='RNA', cell_names='layer3')
  
  celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=sc_ref, sc_assay = 'RNA', 
                                       reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                       dist_thresh=0.55, top_spot=, spot_n=5, repel_r=5, repel_iter=20, keep_model=T)$celltrek
  
  #CellTrek::celltrek_vis(celltrek@meta.data %>% dplyr::select(coord_x, coord_y, layer3:id_new),
  #                       celltrek@images$slice1@image, celltrek@images$slice1@scale.factors$lowres)
  
  write.csv(data.frame(id = celltrek$id_raw, x = celltrek$coord_x, y = celltrek$coord_y), paste("~/data/idCHD/celltrek/", substr(file, 23, 32), "_celltrek5_5.csv", sep = ""), row.names = F)
  

  rm(traint)
  rm(celltrek)
  gc()
}

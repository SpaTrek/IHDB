library(anndata)
library(Matrix)
library(Seurat)

g <- read_h5ad("C:/Users/woloo/Desktop/ST/mouse/raw_h5ad_with_figure/GSM5355663_Sham.h5ad")
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

seu <- readRDS("C:/Users/woloo/Desktop/test.rds")

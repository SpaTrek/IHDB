jaccard <- function (a, b) {
  intersection = length ( intersect (a,b))
  union = length (a) + length (b) - intersection
  return (intersection/union)
}
path_to_robust_programs <- "mouse_spatial_robust_programs.csv"
nmf_programs <- read.csv(path_to_robust_programs, header = T)

mat <- matrix(nrow = length(colnames(nmf_programs)), ncol = length(colnames(nmf_programs)))
for (i in 1:length(colnames(nmf_programs))) {
  for (j in 1:length(colnames(nmf_programs))) {
    mat[i,j] <- jaccard(nmf_programs[[colnames(nmf_programs)[i]]], nmf_programs[[colnames(nmf_programs)[j]]])
  }
}
pheatmap::pheatmap(mat, color = (colorRampPalette(brewer.pal(9,"GnBu"))(100)), clustering_method = "ward.D", cutree_rows = 10, cutree_cols = 10, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", show_rownames = T, fontsize_row = 2)
a <- pheatmap(mat, color = (colorRampPalette(brewer.pal(9,"GnBu"))(100)), clustering_method = "ward.D", cutree_rows = 10, cutree_cols = 10, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
a <- pheatmap::pheatmap(mat, color = (colorRampPalette(brewer.pal(9,"GnBu"))(100)), clustering_method = "ward.D", cutree_rows = 10, cutree_cols = 10, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
View(a)
gc()

order <- a$tree_row$order
samples_ordered <- colnames(nmf_programs)[order]

cluster <- list()
seq <- c(6,3,5,1,7,8,2,9,10)
loc <- 0
for (i in 1:9) {
  cluster[[as.character(i)]] <- samples_ordered[(loc+1):(loc+table(cutree(a$tree_row, k=10))[seq[i]])]
  loc <- loc + table(cutree(a$tree_row, k=10))[seq[i]]
}
cluster

cluster <- list("1"=c('GSM8305655_rank4_9_nruns10.RDS.30.6','GSM8305659_rank4_9_nruns10.RDS.30.2','GSM8305660_rank4_9_nruns10.RDS.30.14',
                      'GSM8305656_rank4_9_nruns10.RDS.30.4','GSM8305653_rank4_9_nruns10.RDS.30.9','GSM8305662_rank4_9_nruns10.RDS.30.6',
                      'GSM8305649_rank4_9_nruns10.RDS.30.7','GSM6613086_rank4_9_nruns10.RDS.30.4','GSM8305643_rank4_9_nruns10.RDS.30.8',
                      'GSM6613081_rank4_9_nruns10.RDS.30.2','GSM6613087_rank4_9_nruns10.RDS.30.20','GSM6613088_rank4_9_nruns10.RDS.30.6',
                      'GSM8305651_rank4_9_nruns10.RDS.30.4','GSM8305635_rank4_9_nruns10.RDS.30.9','GSM8305647_rank4_9_nruns10.RDS.30.12',
                      'GSM8305648_rank4_9_nruns10.RDS.30.2','GSM8305652_rank4_9_nruns10.RDS.30.2','GSM5943195_rank4_9_nruns10.RDS.30.9',
                      'GSM8305640_rank4_9_nruns10.RDS.30.7','GSM8305641_rank4_9_nruns10.RDS.30.6','GSM5943191_rank4_9_nruns10.RDS.30.4',
                      'GSM8305650_rank4_9_nruns10.RDS.30.20','GSM5943198_rank4_9_nruns10.RDS.30.4','GSM8305629_rank4_9_nruns10.RDS.30.14',
                      'GSM5943189_rank4_9_nruns10.RDS.30.3','GSM5943190_rank4_9_nruns10.RDS.30.11','GSM5355666_rank4_9_nruns10.RDS.30.26',
                      'GSM6613089_rank4_9_nruns10.RDS.30.4','GSM8305646_rank4_9_nruns10.RDS.30.12','GSM8305633_rank4_9_nruns10.RDS.30.8',
                      'GSM8305634_rank4_9_nruns10.RDS.30.2','GSM5355666_rank4_9_nruns10.RDS.30.4','GSM6613085_rank4_9_nruns10.RDS.30.22',
                      'GSM6613084_rank4_9_nruns10.RDS.30.11','GSM8305630_rank4_9_nruns10.RDS.30.11','GSM6613086_rank4_9_nruns10.RDS.30.25',
                      'GSM8305661_rank4_9_nruns10.RDS.30.3','GSM8305655_rank4_9_nruns10.RDS.30.14','GSM8305659_rank4_9_nruns10.RDS.30.5',
                      'GSM8305662_rank4_9_nruns10.RDS.30.30','GSM8305657_rank4_9_nruns10.RDS.30.2','GSM8305656_rank4_9_nruns10.RDS.30.21',
                      'GSM8305660_rank4_9_nruns10.RDS.30.13'),
                "2"=c('GSM6613084_rank4_9_nruns10.RDS.30.16','GSM8305630_rank4_9_nruns10.RDS.30.16','GSM6613080_rank4_9_nruns10.RDS.30.11',
                      'GSM5355666_rank4_9_nruns10.RDS.30.5','GSM6613081_rank4_9_nruns10.RDS.30.18','GSM8305641_rank4_9_nruns10.RDS.30.18',
                      'GSM6613086_rank4_9_nruns10.RDS.30.12','GSM8305650_rank4_9_nruns10.RDS.30.14','GSM6613088_rank4_9_nruns10.RDS.30.26',
                      'GSM8305659_rank4_9_nruns10.RDS.30.14','GSM8305633_rank4_9_nruns10.RDS.30.22','GSM8305660_rank4_9_nruns10.RDS.30.11',
                      'GSM8305657_rank4_9_nruns10.RDS.30.16','GSM8305661_rank4_9_nruns10.RDS.30.14','GSM8305655_rank4_9_nruns10.RDS.30.10',
                      'GSM8305656_rank4_9_nruns10.RDS.30.15','GSM8305648_rank4_9_nruns10.RDS.30.20','GSM8305659_rank4_9_nruns10.RDS.30.23',
                      'GSM8305661_rank4_9_nruns10.RDS.30.16','GSM5355663_rank4_9_nruns10.RDS.30.8','GSM5355668_rank4_9_nruns10.RDS.30.16',
                      'GSM5943194_rank4_9_nruns10.RDS.30.4','GSM6613084_rank4_9_nruns10.RDS.30.6','GSM8305630_rank4_9_nruns10.RDS.30.6',
                      'GSM8305650_rank4_9_nruns10.RDS.30.19','GSM5355666_rank4_9_nruns10.RDS.30.18','GSM8305651_rank4_9_nruns10.RDS.30.15',
                      'GSM8305633_rank4_9_nruns10.RDS.30.16','GSM8305634_rank4_9_nruns10.RDS.30.14','GSM6613089_rank4_9_nruns10.RDS.30.15',
                      'GSM8305639_rank4_9_nruns10.RDS.30.20','GSM8305644_rank4_9_nruns10.RDS.30.5','GSM6613078_rank4_9_nruns10.RDS.30.14',
                      'GSM8305642_rank4_9_nruns10.RDS.30.6','GSM6613087_rank4_9_nruns10.RDS.30.2','GSM8305635_rank4_9_nruns10.RDS.30.13',
                      'GSM5943192_rank4_9_nruns10.RDS.30.7','GSM5943193_rank4_9_nruns10.RDS.30.9','GSM6613081_rank4_9_nruns10.RDS.30.4',
                      'GSM6613085_rank4_9_nruns10.RDS.30.29','GSM6613080_rank4_9_nruns10.RDS.30.2','GSM8305634_rank4_9_nruns10.RDS.30.5',
                      'GSM8305645_rank4_9_nruns10.RDS.30.5','GSM8305646_rank4_9_nruns10.RDS.30.5','GSM5943197_rank4_9_nruns10.RDS.30.9',
                      'GSM6613079_rank4_9_nruns10.RDS.30.9','GSM6613083_rank4_9_nruns10.RDS.30.11','GSM6613077_rank4_9_nruns10.RDS.30.11',
                      'GSM6613082_rank4_9_nruns10.RDS.30.12'),
                "3"=c('GSM6613088_rank4_9_nruns10.RDS.30.3','GSM6613089_rank4_9_nruns10.RDS.30.28','GSM8305633_rank4_9_nruns10.RDS.30.4',
                      'GSM8305648_rank4_9_nruns10.RDS.30.18','GSM6613081_rank4_9_nruns10.RDS.30.6','GSM5355666_rank4_9_nruns10.RDS.30.9',
                      'GSM8305647_rank4_9_nruns10.RDS.30.5','GSM5355666_rank4_9_nruns10.RDS.30.2','GSM6613087_rank4_9_nruns10.RDS.30.23',
                      'GSM8305649_rank4_9_nruns10.RDS.30.29','GSM8305653_rank4_9_nruns10.RDS.30.2','GSM6613078_rank4_9_nruns10.RDS.30.5',
                      'GSM8305646_rank4_9_nruns10.RDS.30.6','GSM8305656_rank4_9_nruns10.RDS.30.6','GSM8305657_rank4_9_nruns10.RDS.30.7',
                      'GSM8305662_rank4_9_nruns10.RDS.30.10','GSM6613086_rank4_9_nruns10.RDS.30.16','GSM8305660_rank4_9_nruns10.RDS.30.6',
                      'GSM8305650_rank4_9_nruns10.RDS.30.4','GSM8305640_rank4_9_nruns10.RDS.30.2','GSM8305655_rank4_9_nruns10.RDS.30.20'),
                "4"=c('GSM8305639_rank4_9_nruns10.RDS.30.10','GSM8305643_rank4_9_nruns10.RDS.30.26','GSM5355668_rank4_9_nruns10.RDS.30.7',
                      'GSM5355666_rank4_9_nruns10.RDS.30.10','GSM5943194_rank4_9_nruns10.RDS.30.13','GSM8305659_rank4_9_nruns10.RDS.30.13',
                      'GSM8305645_rank4_9_nruns10.RDS.30.8','GSM8305662_rank4_9_nruns10.RDS.30.9','GSM5943190_rank4_9_nruns10.RDS.30.10',
                      'GSM5943189_rank4_9_nruns10.RDS.30.5','GSM5943198_rank4_9_nruns10.RDS.30.8','GSM8305648_rank4_9_nruns10.RDS.30.21',
                      'GSM8305635_rank4_9_nruns10.RDS.30.10','GSM6613080_rank4_9_nruns10.RDS.30.12','GSM8305633_rank4_9_nruns10.RDS.30.21',
                      'GSM8305629_rank4_9_nruns10.RDS.30.8','GSM5355663_rank4_9_nruns10.RDS.30.3','GSM8305660_rank4_9_nruns10.RDS.30.7',
                      'GSM8305649_rank4_9_nruns10.RDS.30.9','GSM8305653_rank4_9_nruns10.RDS.30.6','GSM8305655_rank4_9_nruns10.RDS.30.13',
                      'GSM5943191_rank4_9_nruns10.RDS.30.7','GSM6613083_rank4_9_nruns10.RDS.30.12','GSM8305644_rank4_9_nruns10.RDS.30.6',
                      'GSM8305634_rank4_9_nruns10.RDS.30.11','GSM8305659_rank4_9_nruns10.RDS.30.6','GSM8305651_rank4_9_nruns10.RDS.30.5',
                      'GSM8305652_rank4_9_nruns10.RDS.30.6','GSM6613084_rank4_9_nruns10.RDS.30.5','GSM8305630_rank4_9_nruns10.RDS.30.5',
                      'GSM8305641_rank4_9_nruns10.RDS.30.4','GSM8305642_rank4_9_nruns10.RDS.30.5','GSM5943196_rank4_9_nruns10.RDS.30.5',
                      'GSM5943199_rank4_9_nruns10.RDS.30.3','GSM5943197_rank4_9_nruns10.RDS.30.2','GSM6613085_rank4_9_nruns10.RDS.30.4',
                      'GSM5355668_rank4_9_nruns10.RDS.30.3','GSM5943195_rank4_9_nruns10.RDS.30.3','GSM8305640_rank4_9_nruns10.RDS.30.24',
                      'GSM8305648_rank4_9_nruns10.RDS.30.10'),
                "5"=c('GSM8305635_rank4_9_nruns10.RDS.30.28','GSM8305647_rank4_9_nruns10.RDS.30.30','GSM6613088_rank4_9_nruns10.RDS.30.2',
                      'GSM8305643_rank4_9_nruns10.RDS.30.12','GSM8305652_rank4_9_nruns10.RDS.30.11','GSM6613089_rank4_9_nruns10.RDS.30.22',
                      'GSM8305639_rank4_9_nruns10.RDS.30.18','GSM6613084_rank4_9_nruns10.RDS.30.2','GSM8305630_rank4_9_nruns10.RDS.30.2',
                      'GSM8305656_rank4_9_nruns10.RDS.30.2','GSM8305633_rank4_9_nruns10.RDS.30.14','GSM6613086_rank4_9_nruns10.RDS.30.2',
                      'GSM8305651_rank4_9_nruns10.RDS.30.2','GSM8305642_rank4_9_nruns10.RDS.30.3','GSM5355666_rank4_9_nruns10.RDS.30.13',
                      'GSM5943190_rank4_9_nruns10.RDS.30.3','GSM8305641_rank4_9_nruns10.RDS.30.10','GSM8305650_rank4_9_nruns10.RDS.30.2'),
                "6"=c('GSM6613088_rank4_9_nruns10.RDS.30.21','GSM5355666_rank4_9_nruns10.RDS.30.25','GSM6613089_rank4_9_nruns10.RDS.30.3',
                      'GSM8305644_rank4_9_nruns10.RDS.30.11','GSM8305647_rank4_9_nruns10.RDS.30.3','GSM8305646_rank4_9_nruns10.RDS.30.3',
                      'GSM8305659_rank4_9_nruns10.RDS.30.11','GSM8305651_rank4_9_nruns10.RDS.30.16','GSM8305653_rank4_9_nruns10.RDS.30.11',
                      'GSM8305661_rank4_9_nruns10.RDS.30.13','GSM8305641_rank4_9_nruns10.RDS.30.8','GSM8305649_rank4_9_nruns10.RDS.30.13',
                      'GSM8305650_rank4_9_nruns10.RDS.30.6','GSM8305635_rank4_9_nruns10.RDS.30.29','GSM8305652_rank4_9_nruns10.RDS.30.16',
                      'GSM6613084_rank4_9_nruns10.RDS.30.3','GSM8305630_rank4_9_nruns10.RDS.30.3','GSM5943189_rank4_9_nruns10.RDS.30.9',
                      'GSM8305629_rank4_9_nruns10.RDS.30.2','GSM6613078_rank4_9_nruns10.RDS.30.4','GSM8305640_rank4_9_nruns10.RDS.30.19',
                      'GSM5943189_rank4_9_nruns10.RDS.30.7','GSM6613089_rank4_9_nruns10.RDS.30.10','GSM8305652_rank4_9_nruns10.RDS.30.4',
                      'GSM8305662_rank4_9_nruns10.RDS.30.8','GSM5943191_rank4_9_nruns10.RDS.30.2','GSM5943190_rank4_9_nruns10.RDS.30.4',
                      'GSM5943198_rank4_9_nruns10.RDS.30.7','GSM8305634_rank4_9_nruns10.RDS.30.4','GSM8305643_rank4_9_nruns10.RDS.30.4',
                      'GSM6613086_rank4_9_nruns10.RDS.30.9','GSM8305655_rank4_9_nruns10.RDS.30.4','GSM8305657_rank4_9_nruns10.RDS.30.4',
                      'GSM8305656_rank4_9_nruns10.RDS.30.3','GSM8305660_rank4_9_nruns10.RDS.30.12'),
                "7"=c('GSM6613078_rank4_9_nruns10.RDS.30.13','GSM8305659_rank4_9_nruns10.RDS.30.9','GSM6613087_rank4_9_nruns10.RDS.30.12',
                      'GSM8305648_rank4_9_nruns10.RDS.30.23','GSM8305655_rank4_9_nruns10.RDS.30.11','GSM6613081_rank4_9_nruns10.RDS.30.23',
                      'GSM8305635_rank4_9_nruns10.RDS.30.12','GSM8305634_rank4_9_nruns10.RDS.30.29','GSM8305643_rank4_9_nruns10.RDS.30.9',
                      'GSM8305640_rank4_9_nruns10.RDS.30.14','GSM8305647_rank4_9_nruns10.RDS.30.18','GSM6613088_rank4_9_nruns10.RDS.30.22',
                      'GSM5355666_rank4_9_nruns10.RDS.30.22','GSM6613089_rank4_9_nruns10.RDS.30.20','GSM8305661_rank4_9_nruns10.RDS.30.15',
                      'GSM8305629_rank4_9_nruns10.RDS.30.9','GSM8305653_rank4_9_nruns10.RDS.30.21','GSM8305646_rank4_9_nruns10.RDS.30.14',
                      'GSM8305656_rank4_9_nruns10.RDS.30.14','GSM8305662_rank4_9_nruns10.RDS.30.17','GSM8305642_rank4_9_nruns10.RDS.30.28',
                      'GSM6613086_rank4_9_nruns10.RDS.30.23','GSM8305651_rank4_9_nruns10.RDS.30.25','GSM8305635_rank4_9_nruns10.RDS.30.5',
                      'GSM8305641_rank4_9_nruns10.RDS.30.23','GSM8305652_rank4_9_nruns10.RDS.30.15','GSM8305646_rank4_9_nruns10.RDS.30.10',
                      'GSM8305657_rank4_9_nruns10.RDS.30.17','GSM8305653_rank4_9_nruns10.RDS.30.27','GSM5943198_rank4_9_nruns10.RDS.30.10',
                      'GSM6613083_rank4_9_nruns10.RDS.30.8','GSM8305629_rank4_9_nruns10.RDS.30.4','GSM5943190_rank4_9_nruns10.RDS.30.7',
                      'GSM5943191_rank4_9_nruns10.RDS.30.18','GSM8305661_rank4_9_nruns10.RDS.30.9','GSM8305633_rank4_9_nruns10.RDS.30.5',
                      'GSM8305655_rank4_9_nruns10.RDS.30.5','GSM8305643_rank4_9_nruns10.RDS.30.21','GSM6613078_rank4_9_nruns10.RDS.30.8',
                      'GSM8305648_rank4_9_nruns10.RDS.30.9','GSM8305649_rank4_9_nruns10.RDS.30.19','GSM8305660_rank4_9_nruns10.RDS.30.5',
                      'GSM5943197_rank4_9_nruns10.RDS.30.4','GSM5943193_rank4_9_nruns10.RDS.30.4','GSM5943194_rank4_9_nruns10.RDS.30.10',
                      'GSM5943192_rank4_9_nruns10.RDS.30.6','GSM5943195_rank4_9_nruns10.RDS.30.11','GSM6613086_rank4_9_nruns10.RDS.30.5',
                      'GSM8305661_rank4_9_nruns10.RDS.30.4','GSM6613084_rank4_9_nruns10.RDS.30.17','GSM8305630_rank4_9_nruns10.RDS.30.17',
                      'GSM6613088_rank4_9_nruns10.RDS.30.8','GSM8305659_rank4_9_nruns10.RDS.30.7','GSM8305662_rank4_9_nruns10.RDS.30.21',
                      'GSM6613077_rank4_9_nruns10.RDS.30.13','GSM8305657_rank4_9_nruns10.RDS.30.12','GSM8305646_rank4_9_nruns10.RDS.30.11',
                      'GSM5355666_rank4_9_nruns10.RDS.30.11','GSM8305647_rank4_9_nruns10.RDS.30.15','GSM6613080_rank4_9_nruns10.RDS.30.22',
                      'GSM8305644_rank4_9_nruns10.RDS.30.14','GSM6613082_rank4_9_nruns10.RDS.30.17','GSM6613083_rank4_9_nruns10.RDS.30.7',
                      'GSM6613078_rank4_9_nruns10.RDS.30.17','GSM6613079_rank4_9_nruns10.RDS.30.15','GSM5355663_rank4_9_nruns10.RDS.30.13',
                      'GSM6613089_rank4_9_nruns10.RDS.30.26','GSM8305657_rank4_9_nruns10.RDS.30.14','GSM6613080_rank4_9_nruns10.RDS.30.23',
                      'GSM8305644_rank4_9_nruns10.RDS.30.18','GSM8305641_rank4_9_nruns10.RDS.30.5','GSM8305650_rank4_9_nruns10.RDS.30.16',
                      'GSM6613078_rank4_9_nruns10.RDS.30.9','GSM6613077_rank4_9_nruns10.RDS.30.26','GSM8305656_rank4_9_nruns10.RDS.30.23',
                      'GSM6613086_rank4_9_nruns10.RDS.30.17','GSM8305655_rank4_9_nruns10.RDS.30.22','GSM8305635_rank4_9_nruns10.RDS.30.20',
                      'GSM8305642_rank4_9_nruns10.RDS.30.11','GSM8305640_rank4_9_nruns10.RDS.30.15','GSM8305651_rank4_9_nruns10.RDS.30.13',
                      'GSM8305645_rank4_9_nruns10.RDS.30.30','GSM8305653_rank4_9_nruns10.RDS.30.29','GSM6613082_rank4_9_nruns10.RDS.30.6',
                      'GSM6613077_rank4_9_nruns10.RDS.30.6','GSM6613083_rank4_9_nruns10.RDS.30.6','GSM8305645_rank4_9_nruns10.RDS.30.3',
                      'GSM8305649_rank4_9_nruns10.RDS.30.16','GSM5943197_rank4_9_nruns10.RDS.30.8','GSM8305655_rank4_9_nruns10.RDS.30.26',
                      'GSM5943196_rank4_9_nruns10.RDS.30.8','GSM5355663_rank4_9_nruns10.RDS.30.7','GSM8305635_rank4_9_nruns10.RDS.30.30',
                      'GSM6613085_rank4_9_nruns10.RDS.30.5','GSM6613081_rank4_9_nruns10.RDS.30.3','GSM8305648_rank4_9_nruns10.RDS.30.12',
                      'GSM6613088_rank4_9_nruns10.RDS.30.14','GSM8305633_rank4_9_nruns10.RDS.30.3','GSM6613087_rank4_9_nruns10.RDS.30.8',
                      'GSM8305629_rank4_9_nruns10.RDS.30.7'),
                "8"=c('GSM8305653_rank4_9_nruns10.RDS.30.8','GSM8305659_rank4_9_nruns10.RDS.30.4','GSM8305660_rank4_9_nruns10.RDS.30.4',
                      'GSM8305662_rank4_9_nruns10.RDS.30.2','GSM8305634_rank4_9_nruns10.RDS.30.7','GSM8305648_rank4_9_nruns10.RDS.30.5',
                      'GSM8305657_rank4_9_nruns10.RDS.30.6','GSM8305645_rank4_9_nruns10.RDS.30.10','GSM8305655_rank4_9_nruns10.RDS.30.2',
                      'GSM6613082_rank4_9_nruns10.RDS.30.7','GSM6613079_rank4_9_nruns10.RDS.30.5','GSM6613083_rank4_9_nruns10.RDS.30.2',
                      'GSM6613085_rank4_9_nruns10.RDS.30.3','GSM6613080_rank4_9_nruns10.RDS.30.17','GSM6613081_rank4_9_nruns10.RDS.30.12',
                      'GSM8305629_rank4_9_nruns10.RDS.30.15','GSM8305649_rank4_9_nruns10.RDS.30.6'),
                "9"=c('GSM6613084_rank4_9_nruns10.RDS.30.18','GSM8305630_rank4_9_nruns10.RDS.30.18','GSM8305645_rank4_9_nruns10.RDS.30.29',
                      'GSM8305643_rank4_9_nruns10.RDS.30.20','GSM8305662_rank4_9_nruns10.RDS.30.25','GSM8305635_rank4_9_nruns10.RDS.30.2',
                      'GSM8305649_rank4_9_nruns10.RDS.30.5','GSM8305640_rank4_9_nruns10.RDS.30.6','GSM8305641_rank4_9_nruns10.RDS.30.7',
                      'GSM8305651_rank4_9_nruns10.RDS.30.6','GSM8305652_rank4_9_nruns10.RDS.30.13','GSM8305642_rank4_9_nruns10.RDS.30.18',
                      'GSM8305650_rank4_9_nruns10.RDS.30.9'))


genes <- list("1"=c(),"2"=c(),"3"=c(),"4"=c(),"5"=c(),"6"=c(),"7"=c(),"8"=c(),"9"=c())
for (i in 1:9) {
  all_genes <- c()
  for (j in 1:length(cluster[[i]])) {
    all_genes <- c(all_genes, nmf_programs[[cluster[[i]][j]]])
  }
  genes[[i]] <- names(table(all_genes)[order(table(all_genes), decreasing = T)])[1:50]
}

jaccard <- function (a, b) {
  intersection = length ( intersect (a,b))
  union = length (a) + length (b) - intersection
  return (intersection/union)
}

mati <- matrix(nrow = 9, ncol = 40)
MP <- read.csv("C:/Users/woloo/华为云盘/Inte-Inte-Inte-Inte-Inte-Inte/idCHD/idCHD投稿/Supplementary Data/Supplementary Data Table 3 mouse_MP_list.csv")
for (i in c(1:9)) {
  for (j in 1:40) {
    mati[i,j] <- jaccard(unlist(genes[i]), MP[,j])
  }
}

MP_used <- c(7, 9, 13,15,17,18,19, 25,27,32,36)
pheatmap::pheatmap(mati[,MP_used], scale = "none", labels_col = MP_used, labels_row = 1:10, color = colorRampPalette(brewer.pal(7, "OrRd"))(100), border_color = "white")

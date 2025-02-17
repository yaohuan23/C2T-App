library(monocle3)
library(magrittr)
library(sva)
setwd("~/Desktop/C2T-seq/Figure3-maternal-RNA/BH1-3/")
Anno <- read.csv("./geneAnno.csv",header = TRUE,row.names = 1)
my_sc_expression <-read.csv("PE150-250-tail.csv",header = TRUE,row.names = 1)
sum <- as.data.frame(colSums(my_sc_expression))
mata_data <- read.csv("ABH13-meta.csv",header = TRUE,row.names = 1,check.names = F)
my_sc_expression <-as.matrix(my_sc_expression)
rownames(mata_data) <- colnames(my_sc_expression)
mata_data$RNAMount <- sum$`colSums(my_sc_expression)`
Myexp <-new_cell_data_set(my_sc_expression,cell_metadata = mata_data )
#Myexp <- preprocess_cds(Myexp, num_dim = 10,norm_method = "log")
Myexp <- preprocess_cds(Myexp, num_dim = 10, norm_method = "size_only")
plot_pc_variance_explained(Myexp)

#threedcds <- reduce_dimension(cds,max_components = 100)
Myexp <- reduce_dimension(Myexp,reduction_method="UMAP")
Myexp <- cluster_cells(Myexp,resolution=0.05)
clus_info<- clusters(Myexp)
Cell_cinfo <- as.numeric(clus_info[row.names(mata_data)])
mata_data <- cbind(mata_data,Cell_cinfo)


plot_cells(Myexp,cell_size = 1.5)
plot_cells(Myexp,color_cells_by = "Sample",cell_size = 1)
plot_cells(Myexp,color_cells_by = "geneType",cell_size = 1)
plot_cells(Myexp,color_cells_by = "TransIn",cell_size = 1)
plot_cells(Myexp,color_cells_by = "Lable",cell_size = 1)
#GATA
plot_cells(Myexp, genes=c("GATA2","DLX3"),show_trajectory_graph = FALSE,norm_method = "log",cell_size = 0.5)
plot_cells(Myexp, genes=c("DLX3"),show_trajectory_graph = FALSE,norm_method = "log",cell_size = 1)

mod_matrix <- model.matrix(~ batch + Lable, data = colData(Myexp))
corrected_data <- ComBat(dat = counts(Myexp), batch = colData(Myexp)$batch, mod = mod_matrix)












# non-Cap data analysis

Anno <- read.csv("./geneAnno.csv",header = TRUE,row.names = 1)
my_sc_expression <-read.csv("Non-Cap-m7G.csv",header = TRUE,row.names = 1)
sum <- as.data.frame(colSums(my_sc_expression))
mata_data <- read.csv("ABH13-meta.csv",header = TRUE,row.names = 1,check.names = F)
my_sc_expression <-as.matrix(my_sc_expression)
rownames(mata_data) <- colnames(my_sc_expression)
mata_data$RNAMount <- sum$`colSums(my_sc_expression)`
Myexp <-new_cell_data_set(my_sc_expression,cell_metadata = mata_data)
Myexp <- preprocess_cds(Myexp, num_dim = 8,norm_method = "log")
#Myexp <- preprocess_cds(Myexp, num_dim = 10, norm_method = "size_only")
plot_pc_variance_explained(Myexp)


Myexp <- reduce_dimension(Myexp,reduction_method="UMAP")
Myexp <- cluster_cells(Myexp,resolution=0.05)
clus_info<- clusters(Myexp)
Cell_cinfo <- as.numeric(clus_info[row.names(mata_data)])
mata_data <- cbind(mata_data,Cell_cinfo)
plot_cells(Myexp,cell_size = 1)
plot_cells(Myexp,color_cells_by = "Sample",cell_size = 1)
plot_cells(Myexp,color_cells_by = "geneType",cell_size = 1)
plot_cells(Myexp,color_cells_by = "TransIn",cell_size = 1)
plot_cells(Myexp,color_cells_by = "Lable",cell_size = 1)
plot_cells(Myexp,color_cells_by = "batch",cell_size = 1)

# Cap data analysis

Anno <- read.csv("./geneAnno.csv",header = TRUE,row.names = 1)
my_sc_expression <-read.csv("Cap-m7G.csv",header = TRUE,row.names = 1)
mata_data <- read.csv("ABH13-meta.csv",header = TRUE,row.names = 1,check.names = F)
my_sc_expression <-as.matrix(my_sc_expression)
rownames(mata_data) <- colnames(my_sc_expression)
Myexp <-new_cell_data_set(my_sc_expression,cell_metadata = mata_data)
#Myexp <- preprocess_cds(Myexp, num_dim = 10,norm_method = "log")
Myexp <- preprocess_cds(Myexp, num_dim = 7, norm_method = "size_only")
plot_pc_variance_explained(Myexp)


Myexp <- reduce_dimension(Myexp,reduction_method="UMAP")
Myexp <- cluster_cells(Myexp,resolution=0.05)
clus_info<- clusters(Myexp)
Cell_cinfo <- as.numeric(clus_info[row.names(mata_data)])
mata_data <- cbind(mata_data,Cell_cinfo)
plot_cells(Myexp,cell_size = 1)
plot_cells(Myexp,color_cells_by = "Sample",cell_size = 1)
plot_cells(Myexp,color_cells_by = "geneType",cell_size = 1)
plot_cells(Myexp,color_cells_by = "TransIn",cell_size = 1)
plot_cells(Myexp,color_cells_by = "Lable",cell_size = 1)
plot_cells(Myexp,color_cells_by = "batch",cell_size = 1)

#####################################

# merged_analysis
my_sc_expression_tail <-read.csv("PE150-250-tail.csv",header = TRUE,row.names = 1)
my_sc_expression_noCap <-read.csv("Non-Cap-m7G.csv",header = TRUE,row.names = 1)
my_sc_expression_Cap <-read.csv("Cap-m7G.csv",header = TRUE,row.names = 1)
my_sc_expression_noCap <-log2(my_sc_expression_noCap+1)
variances <- apply(my_sc_expression_tail, 1, var)
nonzero_columns <- variances != 0
pca_tail <- prcomp(t(my_sc_expression_tail[nonzero_columns,]), scale. = TRUE)
reduced_matrix_tail <- as.matrix(pca_tail$x[, 1:2])
variances <- apply(my_sc_expression_noCap, 1, var)
nonzero_columns <- variances != 0
pca_noCap <- prcomp(t(my_sc_expression_noCap[nonzero_columns,]), scale. = TRUE)
reduced_matrix_noCap <- as.matrix(pca_noCap$x[, 2:3])
variances <- apply(my_sc_expression_Cap, 1, var)
nonzero_columns <- variances != 0
pca_Cap <- prcomp(t(my_sc_expression_Cap[nonzero_columns,]), scale. = TRUE)
reduced_matrix_Cap <- as.matrix(pca_Cap$x[, 2:3])
conbid_matrix <- cbind(reduced_matrix_noCap,reduced_matrix_Cap,reduced_matrix_tail)

umap_result <- umap(conbid_matrix)
reduced_coordinates <- umap_result$layout
df <- data.frame(ID = 1:nrow(reduced_coordinates), 
                 X = reduced_coordinates[, 1],
                 Y = reduced_coordinates[, 2])
df$color <- mata_data$Lable
ggplot(df, aes(x = X, y = Y,color = df$color)) +
  geom_point() +
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
  ggtitle("UMAP Visualization")



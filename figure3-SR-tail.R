library(monocle3)
library(magrittr)
library()
setwd("~/Desktop/C2T-seq/Figure3-maternal-RNA/BH1-3/SR-tail/")
log2_NonCap <- read.table("./log2-NonCap_GV_MII_IN.tsv",row.names = 1)
log2_Cap <- read.table("./log2-CPM-Cap_GV_MII_IN.tsv",row.names = 1)
log2_tail <- read.table("./log2-SRTail-GV-MII-MII_IN.tsv",row.names = 1)
log2_Tanslation <- read.table("./log2-RPKM-merge-xieweiTrans.tsv",row.names = 1)
log2_Tanslation <- log2_Tanslation[,1:2]

row_sums <- rowSums(log2_NonCap)
log2_NonCap <- log2_NonCap[row_sums > 1, ]
log2_NonCap <- log2_NonCap[intersect(rownames(log2_Tanslation),rownames(log2_NonCap)), ]

a<- rownames(log2_NonCap)
log2_Cap <-log2_Cap[a,]
log2_tail <-log2_tail[a,]
log2_Tanslation <- log2_Tanslation[a,]

plot(log2_NonCap$V2,log2_Cap$V2)
plot(log2_NonCap$V2,log2_tail$V2)
plot(log2_Cap$V2,log2_tail$V2)

GV_data <- cbind(log2_NonCap$V2,log2_Cap$V2,log2_tail$V2,log2_Tanslation$V2)
colnames(GV_data) <- c("body","Cap","tail","translation")
row.names(GV_data) <- row.names(log2_NonCap)
GV_mata_data <- as.matrix(colnames(GV_data))
rownames(GV_mata_data) <-GV_mata_data

MII_data <- cbind(log2_NonCap$V3,log2_Cap$V3,log2_tail$V3,log2_Tanslation$V2)
colnames(MII_data) <- c("MII_body","MII_Cap","MII_tail","MII_translation")
row.names(MII_data) <- row.names(log2_NonCap)
MII_mata_data <- as.matrix(colnames(MII_data))
rownames(MII_mata_data) <-MII_mata_data


geneData <- cbind(GV_data,MII_data)
geneData <- t(geneData)
gene_mata_data <- as.matrix(colnames(geneData))
rownames(gene_mata_data) <- gene_mata_data
geneData_cds <-new_cell_data_set(geneData,cell_metada= gene_mata_data)
geneData_cds <- geneData_cds[,Matrix::colSums(exprs(geneData_cds)) != 0] 
geneData_cds <- preprocess_cds(geneData_cds, num_dim = 2, norm_method = "size_only")
plot_pc_variance_explained(geneData_cds)

#threedcds <- reduce_dimension(cds,max_components = 100)
geneData_cds <- reduce_dimension(geneData_cds,reduction_method="PCA")
geneData_cds <- cluster_cells(geneData_cds,resolution=2)
plot_cells(geneData_cds,cell_size = 1.5)








library(FactoMineR)
pca_result <- prcomp(t(geneData), scale = TRUE)
reduced_data <- pca_result$x[, 1:2]
plot(reduced_data,cex = 0.5)
library(umap)
umap_result <- umap::umap(pca_result$x, n_components = 5)
plot(umap_result$layout[,1], umap_result$layout[,2], asp = 1, xlab = "UMAP 1", ylab = "UMAP 2",cex=0.1,pch=20)
library(dbscan)

dbscan_result <- dbscan(umap_result$layout, eps = 1, minPts = 10)
plot(umap_result$layout, col = dbscan_result$cluster + 1L,cex=0.1,pch=20)
umap_result$class <- dbscan_result$cluster
my_result <- cbind(umap_result$layout,umap_result$class)
write.csv(my_result,"GV-MII-gene-classes.csv")
my_result <-as.data.frame(my_result)


log2_NonCap$class <- my_result$V3
log2_Cap$class <- my_result$V3
log2_tail$class <- my_result$V3
sorted_indices <- order(log2_NonCap[, 4])
log2_NonCap <- log2_NonCap[sorted_indices, ]
sorted_indices <- order(log2_Cap[, 4])
log2_Cap <- log2_Cap[sorted_indices, ]
sorted_indices <- order(log2_tail[, 4])
log2_tail <- log2_tail[sorted_indices, ]

class5_log2_NonCap <-subset(log2_NonCap,log2_NonCap$class==5)
class5_log2_Cap <-subset(log2_Cap,log2_Cap$class==5)
class5_log2_tail <-subset(log2_tail,log2_tail$class==5)
plot(class5_log2_tail$V2,class5_log2_tail$V3)
plot(class5_log2_Cap$V2,class5_log2_Cap$V3)
plot(class5_log2_Cap$V2,class5_log2_Cap$V3)

plot(log2_NonCap$V2,log2_Cap$V2,col=log2_NonCap$class + 1L ,cex=0.5,pch=20)
plot(log2_NonCap$V2,log2_tail$V2,col=log2_NonCap$class + 1L ,cex=0.5,pch=20)

plot(log2_NonCap$V3,log2_Cap$V3,col=log2_NonCap$class + 1L ,cex=0.5,pch=20)
plot(log2_NonCap$V3,log2_tail$V3,col=log2_NonCap$class + 1L ,cex=0.5,pch=20)

plot(log2_tail$V2,log2_tail$V4,col=log2_tail$class +1L,cex=0.5,pch=20)
plot(log2_Cap$V2,log2_Cap$V3,col=log2_tail$class +1L,cex=0.5,pch=20)
plot(log2_Cap$V2,log2_Cap$V4,col=log2_tail$class +1L,cex=0.5,pch=20)
plot(log2_NonCap$V2,log2_NonCap$V3,col=log2_tail$class +1L,cex=0.5,pch=20)

plot(log2_tail$V2,log2_tail$V3,col=log2_tail$class+ 1L  ,cex=0.5,pch=20)
plot(log2_tail$V2,log2_tail$V4,col=log2_tail$class+ 1L  ,cex=0.5,pch=20)
plot(log2_Cap$V2,log2_Cap$V3,col=log2_tail$class+ 1L ,cex=0.5,pch=20)
plot(log2_NonCap$V2,log2_NonCap$V3,col=log2_tail$class+ 1L ,cex=0.5,pch=20)

plot(2*log2_NonCap$V2 - log2_Cap$V2 -log2_tail$V2,log2_NonCap$V3 -log2_NonCap$V2,cex=0.2)



########GV_along

if (!require("FactoMineR")) install.packages("FactoMineR")
library(FactoMineR)
pca_result <- prcomp(t(GV_data), scale = TRUE)
reduced_data <- pca_result$x[, 1:2]
plot(reduced_data,cex = 0.5)
set.seed(123) 
kmeans_result <- kmeans(pca_result$x, centers=5)
print(kmeans_result)
plot(reduced_data, col = kmeans_result$cluster,cex=0.5,pch=20)
write.csv(kmeans_result$cluster,"GV_genetypes.csv")
GV_cluster <- as.data.frame(kmeans_result$cluster)
GV_data <- cbind(log2_NonCap$V2,log2_Cap$V2,log2_tail$V2)
GV_data <-as.data.frame(GV_data)
GV_data$class <- GV_cluster$`kmeans_result$cluster`
colnames(GV_data) <- c("body","Cap","tail","class")
rownames(GV_data) <- rownames(log2_Cap)
plot(GV_data$body,GV_data$tail, col = GV_data$class,cex=0.5,pch=20)
plot(GV_data$body,GV_data$Cap, col = GV_data$class,cex=0.5,pch=20)
plot(GV_data$tail,GV_data$Cap, col = GV_data$class,cex=0.5,pch=20)

log2_NonCap$class <- GV_cluster$`kmeans_result$cluster`
log2_Cap$class <- GV_cluster$`kmeans_result$cluster`
log2_tail$class <- GV_cluster$`kmeans_result$cluster`
plot(log2_NonCap$V2,log2_NonCap$V3,col=log2_NonCap$class,cex=0.5,pch=20)
plot(log2_tail$V2,log2_tail$V3,col=log2_tail$class,cex=0.5,pch=20)


GV_data <- cbind(log2_NonCap$V2,log2_Cap$V2,log2_tail$V2)
colnames(GV_data) <- c("body","Cap","tail")
row.names(GV_data) <- row.names(log2_NonCap)
GV_mata_data <- as.matrix(colnames(GV_data))
rownames(GV_mata_data) <-GV_mata_data
if (!require("FactoMineR")) install.packages("FactoMineR")
library(FactoMineR)
pca_result <- prcomp(t(GV_data), scale = TRUE)
reduced_data <- pca_result$x[, 1:2]
plot(reduced_data,cex = 0.5)
set.seed(123) 
kmeans_result <- kmeans(pca_result$x, centers=5)
print(kmeans_result)
plot(reduced_data, col = kmeans_result$cluster,cex=0.5,pch=20)
write.csv(kmeans_result$cluster,"GV_genetypes.csv")
GV_cluster <- as.data.frame(kmeans_result$cluster)

library(monocle3)
library(magrittr)
library(dplyr)
my_sc_expression_A <- read.csv("PolyA-counts.csv",header = TRUE,row.names = 1)
my_sc_expression_A <- as.matrix(my_sc_expression_A)

my_sc_expression_Cap <- read.csv("Capcounts.csv",header = TRUE,row.names = 1)
my_sc_expression_Cap <- as.matrix(my_sc_expression_Cap)

my_sc_expression_GeneBody <- read.csv("body-counts.csv",header = TRUE,row.names = 1)
my_sc_expression_GeneBody <- as.matrix(my_sc_expression_GeneBody)


mata_data <- read.csv("mata_data.csv",header = TRUE,row.names = 1,check.names = F)
Myexp_A <- new_cell_data_set(my_sc_expression_A,cell_metadata = mata_data )
Myexp_Cap <- new_cell_data_set(my_sc_expression_Cap,cell_metadata = mata_data )
Myexp_GeneBody <- new_cell_data_set(my_sc_expression_GeneBody,cell_metadata = mata_data )


#Myexp <- preprocess_cds(Myexp, num_dim = 10,norm_method = "log")
Myexp_A <- preprocess_cds(Myexp_A, num_dim = 8, norm_method = "size_only")
plot_pc_variance_explained(Myexp_A)
Myexp_A <- reduce_dimension(Myexp_A,reduction_method="PCA")
plot_cells(Myexp_A,cell_size = 1,color_cells_by = "Type",reduction_method ="PCA")
plot_cells(Myexp_A,cell_size = 1,color_cells_by = "Defined",reduction_method ="PCA")



Myexp_Cap <- preprocess_cds(Myexp_Cap, num_dim = 8, norm_method = "size_only")
plot_pc_variance_explained(Myexp_Cap)
Myexp_Cap <- reduce_dimension(Myexp_Cap,reduction_method="PCA")
plot_cells(Myexp_Cap,cell_size = 1,color_cells_by = "Type",reduction_method ="PCA")


Myexp_GeneBody <- preprocess_cds(Myexp_GeneBody, num_dim = 8, norm_method = "size_only")
plot_pc_variance_explained(Myexp_GeneBody)
Myexp_GeneBody <- reduce_dimension(Myexp_GeneBody,reduction_method="PCA")
plot_cells(Myexp_GeneBody,cell_size = 1,color_cells_by = "Type",reduction_method ="PCA")

PCA_A<-as.matrix(reducedDims(Myexp_A)[["PCA"]])
PCA_Cap<-as.matrix(reducedDims(Myexp_Cap)[["PCA"]])
PCA_GeneBody<-as.matrix(reducedDims(Myexp_GeneBody)[["PCA"]])

library(e1071)
PCA_A <- as.data.frame(PCA_A[,1:2])
PCA_A$class <- mata_data$Type
#plot_cells_3d(b[,1],b[,2],b[,3])
spline_PCA_A <- subset(PCA_A,PCA_A$class != "A")
spline_PCA_A <- subset(spline_PCA_A,spline_PCA_A$class != "GC")
spline_matedata_1 <- as.data.frame(spline_PCA_A$class)
rownames(spline_matedata_1) <- c(rownames(spline_PCA_A))
colnames(spline_matedata_1) <- "class_info"
spline_matedata_1 <- spline_PCA_A %>% mutate_at(vars(class),list(~ifelse(. == "GV",0,1)))
names(spline_matedata_1) <- c("x","y","class_info")
svm_data <- spline_matedata_1

write.csv(svm_data,"./A-svm.csv")
svm_data$class_info <- as.factor(svm_data$class_info)
svm_data <- read.csv("./A-svm.csv")
###python###
#import pandas as pd
#from sklearn.svm import SVC
#import numpy as np

# 加载数据
#file_path = '/path/to/your/file.csv'  # 替换为您的文件路径
#data = pd.read_csv(file_path)

# 去除非特征列
#data_cleaned = data.drop(columns=['Unnamed: 0'])

# 分割数据为特征和目标变量
#X = data_cleaned[['x', 'y']]
#y = data_cleaned['class_info']

# 训练线性 SVM 模型
#svm_model_linear = SVC(kernel='linear')
#svm_model_linear.fit(X, y)

# 提取决策边界的参数
#w = svm_model_linear.coef_[0]
#b = svm_model_linear.intercept_[0]

# 显示权重和截距
#print("Weights (Coefficients):", w)
#print("Intercept:", b)

# 计算决策边界的直线方程
#slope = -w[0] / w[1]
#intercept = -b / w[1]
#print("Decision Boundary Line: y = {:.2f} * x + {:.2f}".format(slope, intercept))
#########################
model <- svm(class_info ~ ., data = svm_data, kernel = "linear")
x <- seq(min(svm_data$x), max(svm_data$x), length = 100)
y <- (-0.00301397 * x + 5.626849252297937) / -0.45596159
svm_linear <- data.frame(x = x, y = y)

ggplot(PCA_A,aes(x=PC1,y=PC2)) +geom_point(aes(color =class),size =2) + geom_line(data = svm_linear, aes(x=x,y=y),color ="black")

PCA_All <- cbind(PCA_A[,1:2],PCA_Cap[,1:2],PCA_GeneBody[,1:2])
write.csv(PCA_All, "./merged_PCA.csv")

library(FactoMineR)
pca_result <- prcomp(PCA_All, scale = TRUE)
reduced_data <- pca_result$x[, 1:2]
reduced_data <- as.data.frame(reduced_data)
reduced_data$Type <- mata_data$Type
unique_types <- unique(reduced_data$Type)
colors <- c("red","green","blue","black")
color_mapping <- setNames(colors, unique_types)
plot(reduced_data$PC1,reduced_data$PC2,col=color_mapping[reduced_data$Type],cex = 1, asp =1,pch=20)

library(umap)
umap_result <- umap::umap(pca_result$x[,1:6], n_components = 2)
plot(umap_result$layout[,1], umap_result$layout[,2], asp = 1, xlab = "UMAP 1", ylab = "UMAP 2",cex=1,pch=20)











library(dbscan)


heatmap(PCA_All)



write.csv(b,"PCA-matrix.csv")


my_sc_expression <- as.matrix(m7G_counts)
#sum <- as.data.frame(colSums(my_sc_expression))
mata_data <- read.csv("mata_data.csv",header = TRUE,row.names = 1,check.names = F)
Myexp <-new_cell_data_set(my_sc_expression,cell_metadata = mata_data )

#Myexp <- preprocess_cds(Myexp, num_dim = 10,norm_method = "log")
Myexp <- preprocess_cds(Myexp, num_dim = 8, norm_method = "size_only")
plot_pc_variance_explained(Myexp)
Myexp <- reduce_dimension(Myexp,reduction_method="PCA")
plot_cells(Myexp,cell_size = 2,reduction_method ="PCA")
plot_cells(Myexp,color_cells_by = "Stage",cell_size = 1.5,reduction_method ="PCA")



a<- read.table("polyA-log2CPM-merged.txt")
rownames(a) <- a$V1
Merge_experssion_matrix <- a[,2:6]
colnames(Merge_experssion_matrix) <- c("GV","MII","MI-1","MI-2","GC")
Merge_experssion_matrix <- subset(Merge_experssion_matrix,rowSums(Merge_experssion_matrix) > 4)
# 假设Merge_experssion_matrix已经按照您之前的R代码加载和准备好
# 加载必要的库
library(ggplot2)
library(pheatmap)

# 检查并处理NA值（这里是一个简单的示例，您可能需要根据您的数据情况调整）
Merge_experssion_matrix[is.na(Merge_experssion_matrix)] <- 0

# 标准化数据（Z分数归一化）
standardized_data <- as.data.frame(scale(Merge_experssion_matrix))

# 创建热图
pheatmap(standardized_data, 
         color = colorRampPalette(c("blue", "white", "red"))(255), 
         border_color = NA)

Oocyte_specific <- read.csv("Oocyte_specific_genes.csv")
rownames(Oocyte_specific) <- Oocyte_specific$Gene
Oocyte_specific <- Oocyte_specific[,2:6]
Oocyte_specific <-as.matrix(Oocyte_specific)
write.table(Oocyte_specific,"./Oocyte_specific_genes.tsv")
heatmap(Oocyte_specific)









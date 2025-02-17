library(monocle3)
library(magrittr)
library(dplyr)
#my_sc_expression_Total <- read.csv("Total-normal-reads-single-cell.csv",header = TRUE,row.names = 1)
#my_sc_expression_Total <- as.matrix(my_sc_expression_Total)
#my_log <- log(my_sc_expression_Total+1)
#my_log <- scale(my_log)
#heatmap(my_log)



my_sc_expression_A <- read.csv("PolyA-counts.csv",header = TRUE,row.names = 1)
my_sc_expression_A <- as.matrix(my_sc_expression_A)
my_sc_expression_A <- scale(my_sc_expression_A)
heatmap(my_sc_expression_A)

my_sc_expression_Cap <- read.csv("Capcounts.csv",header = TRUE,row.names = 1)
my_sc_expression_Cap <- as.matrix(my_sc_expression_Cap)
my_sc_expression_Cap <- scale(my_sc_expression_Cap)
heatmap(my_sc_expression_Cap)


my_sc_expression_GeneBody <- read.csv("body-counts.csv",header = TRUE,row.names = 1)
my_sc_expression_GeneBody <- as.matrix(my_sc_expression_GeneBody)
my_sc_expression_GeneBody<- log(my_sc_expression_GeneBody+1)
my_sc_expression_GeneBody <- scale(my_sc_expression_GeneBody)
heatmap(my_sc_expression_GeneBody)

m6A_my_sc_expression_Cap <- read.csv("Cap-m6A-counts.csv",header = TRUE,row.names = 1)
m6A_my_sc_expression_Cap <- as.matrix(m6A_my_sc_expression_Cap)
m6A_my_sc_expression_Cap <- scale(m6A_my_sc_expression_Cap)
heatmap(m6A_my_sc_expression_Cap)

m6A_my_sc_expression_GeneBody <- read.csv("body-m6A-counts.csv",header = TRUE, row.names = 1)
m6A_my_sc_expression_GeneBody <- as.matrix(m6A_my_sc_expression_GeneBody)
m6A_my_sc_expression_GeneBody <- log(m6A_my_sc_expression_GeneBody+1)
m6A_my_sc_expression_GeneBody <- scale(m6A_my_sc_expression_GeneBody)
heatmap(m6A_my_sc_expression_GeneBody)

PCA_A <- t(my_sc_expression_A)
constant_columns <- apply(PCA_A, 2, function(col) length(unique(col)) == 1)
PCA_A_clean <- PCA_A[, !constant_columns]
PCA_A <- prcomp(PCA_A_clean, center = TRUE, scale. = TRUE)
PCA_A_matrix <- PCA_A$x

PCA_Cap <- t(my_sc_expression_Cap)
constant_columns <- apply(PCA_Cap, 2, function(col) length(unique(col)) == 1)
PCA_Cap_clean <- PCA_Cap[, !constant_columns]
PCA_Cap <- prcomp(PCA_Cap_clean, center = TRUE, scale. = TRUE)
PCA_Cap_matrix <- PCA_Cap$x

PCA_body <- t(my_sc_expression_GeneBody)
constant_columns <- apply(PCA_body, 2, function(col) length(unique(col)) == 1)
PCA_body_clean <- PCA_body[, !constant_columns]
PCA_body <- prcomp(PCA_body_clean, center = TRUE, scale. = TRUE)
PCA_body_matrix <- PCA_body$x

PCA_m6A_cap <- t(m6A_my_sc_expression_Cap)
constant_columns <- apply(PCA_m6A_cap, 2, function(col) length(unique(col)) == 1)
PCA_m6A_Cap_clean <- PCA_m6A_cap[, !constant_columns]
PCA_m6A_cap <- prcomp(PCA_m6A_Cap_clean, center = TRUE, scale. = TRUE)
PCA_m6A_Cap_matrix <- PCA_m6A_cap$x


PCA_m6A_body <- t(m6A_my_sc_expression_GeneBody)
constant_columns <- apply(PCA_m6A_body, 2, function(col) length(unique(col)) == 1)
PCA_m6A_body_clean <- PCA_m6A_body[, !constant_columns]
PCA_m6A_body <- prcomp(PCA_m6A_body_clean, center = TRUE, scale. = TRUE)
PCA_m6A_body_matrix <- PCA_m6A_body$x

PCA_A_len <- t(A_tail_150)
constant_columns <- apply(PCA_A_len, 2, function(col) length(unique(col)) == 1)
PCA_A_len_clean <- PCA_A_len[, !constant_columns]
PCA_A_len <- prcomp(PCA_A_len_clean, center = TRUE, scale. = TRUE)
PCA_A_len_matrix <- PCA_A_len$x


A_tail_150 <- read.table("PE150-A-result.txt")
A_tail_150 <- as.matrix(A_tail_150)
heatmap(A_tail_150)


A_tail_250 <- read.table("PE250-A-length.txt")
A_tail_250 <- as.matrix(A_tail_250)
heatmap(A_tail_250)


merged_tail <- read.table("merged-150-250-length.txt",row.names = 1)
colnames(merged_tail) <- c("PE150","PE250")
plot(merged_tail$PE150,merged_tail$PE250,pch = 19, cex = 0.5, col = "darkblue",xlim=c(0,150),ylim=c(0,150))
abline(a = 0, b = 1, col = "red") 
cor(merged_tail$PE150,merged_tail$PE250)
abline(a = -10, b = 1, col = "red") 

APA_150 <- rownames(subset(A_tail_150,A_tail_150[,"ACATCG"]>0))
APA_250 <- rownames(subset(A_tail_250,A_tail_250[,"ACATCG"]>0))
common_APA <- intersect(APA_250,APA_150)
APA_150 <- A_tail_150[common_APA,]
APA_250 <- A_tail_250[common_APA,]

plot(APA_250[,"ACATCG"],APA_150[,"ACATCG"],pch = 19, cex = 0.5, col = "darkblue",xlim=c(0,150),ylim=c(0,150))
abline(a = 0, b = 1, col = "red") 
t.test(APA_250[,"ACATCG"],APA_150[,"ACATCG"],paired = TRUE)





APA_150 <- subset(A_tail_150,A_tail_150[,"ACATCG"]>0)
APA_150 <- subset(APA_150,APA_150[,"GTAGCC"]>0)
plot(APA_150[,"ACATCG"],APA_150[,"GTAGCC"],pch = 19, cex = 0.5, col = "darkblue",xlim=c(0,150),ylim=c(0,150))
abline(a = 0, b = 1, col = "red") 

cor(APA_150[,"ACATCG"],APA_150[,"GTAGCC"])




multiOmic_PCA <- cbind(PCA_A_matrix[,1:2],PCA_A_len_matrix[,1:2],PCA_body_matrix[,1:2],PCA_Cap_matrix[,1:2],PCA_m6A_body_matrix[,1:2],PCA_m6A_Cap_matrix[,1:2])

multiOmic_PCA_2 <- prcomp(multiOmic_PCA, center = TRUE, scale. = TRUE)
multiOmic_PCA_2 <- multiOmic_PCA_2$x
multiOmic_PCA_2 <- multiOmic_PCA_2[,1:2]

write.csv(multiOmic_PCA_2,"All-merged-PCA2.csv")
plot(multiOmic_PCA_2[,1],multiOmic_PCA_2[,2])

write.csv(multiOmic_PCA,"PCA-merge.csv")
multiOmic_PCA <- t(multiOmic_PCA)
heatmap(multiOmic_PCA)

set.seed(1)
umap_result <- umap(multiOmic_PCA, n_components = 3)


my_data <- read.csv("./All-merged-PCA2.csv")

ggplot(data = my_data, aes(x = PC1, y = PC2, color = Type,size = 1)) +
  geom_point()










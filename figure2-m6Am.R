# this script was used to calculate the repeat performence of m7G calling on public data and self data
library(monocle3)
library(magrittr)
setwd("~/Desktop/C2T-seq/Figure2-m6Am")
#non_cap reads
merged_no_Cap <- read.csv("qua_matrix.csv",row.names = 1)
filter_matrix_by_column_sum <- function(matrix, target_sum) {
  column_sums <- colSums(matrix)
  target_columns <- which(column_sums > target_sum)
  filtered_matrix <- matrix[, target_columns]
  return(filtered_matrix)
}
mata_data <- read.csv("GVcKOWT-noCap-meta.csv",header = TRUE,row.names = 1,check.names = F)
expression_data <- merged_no_Cap[,row.names(mata_data)]
expression_data <-as.matrix(expression_data)
Myexp <-new_cell_data_set(expression_data,cell_metadata = mata_data)
Myexp <- preprocess_cds(Myexp, num_dim = 18, norm_method = "log")
plot_pc_variance_explained(Myexp)
Myexp <- reduce_dimension(Myexp,reduction_method="UMAP")
Myexp <- cluster_cells(Myexp,resolution=0.05)
clus_info<- clusters(Myexp)
Cell_cinfo <- as.numeric(clus_info[row.names(mata_data)])
mata_data <- cbind(mata_data,Cell_cinfo)
plot_cells(Myexp,cell_size = 1.5)
plot_cells(Myexp,color_cells_by = "Genome_type",cell_size = 1)

# m7G reads
merged_Cap <- read.csv("qua-GVcKOWT-m7G-matrix.csv",row.names = 1)
expression_data <- merged_Cap[,row.names(mata_data)]
expression_data <-as.matrix(expression_data)
Myexp <-new_cell_data_set(expression_data,cell_metadata = mata_data)
Myexp <- preprocess_cds(Myexp, num_dim = 12, norm_method = "size_only")
plot_pc_variance_explained(Myexp)
Myexp <- reduce_dimension(Myexp,reduction_method="UMAP")
Myexp <- cluster_cells(Myexp,resolution=0.05)
clus_info<- clusters(Myexp)
Cell_cinfo <- as.numeric(clus_info[row.names(mata_data)])
mata_data <- cbind(mata_data,Cell_cinfo)
plot_cells(Myexp,cell_size = 1.5)
plot_cells(Myexp,color_cells_by = "Genome_type",cell_size = 1)


#m6Am reads
merged_m6A <- read.csv("GVcKOWT-m6Am.csv",row.names = 1)
expression_data <- merged_m6A[,row.names(mata_data)]
write.csv(expression_data,"./GVcKO-m6Am.csv")
expression_data <-as.matrix(expression_data)
Myexp <-new_cell_data_set(expression_data,cell_metadata = mata_data)
Myexp <- preprocess_cds(Myexp, num_dim = 10, norm_method = "size_only")
plot_pc_variance_explained(Myexp)
Myexp <- reduce_dimension(Myexp,reduction_method="UMAP")
Myexp <- cluster_cells(Myexp,resolution=0.05)
clus_info<- clusters(Myexp)
Cell_cinfo <- as.numeric(clus_info[row.names(mata_data)])
mata_data <- cbind(mata_data,Cell_cinfo)
plot_cells(Myexp,cell_size = 1.5)
plot_cells(Myexp,color_cells_by = "Genome_type",cell_size = 1)

# Find the m6Am and Cap genes
my_expression <- read.table("./noCap-Cap-m6Am.RPKM.txt")
#caculation of the correlation of repeats
my_expression <-as.matrix(my_expression)
a<-my_expression
a$V2=log2(a$V2+1)
a$V3=log2(a$V3+1)
a$V4=log2(a$V4+1)
a$diff <- a$V4-a$V2
a$diff <- ifelse(a$diff > -2.4 , 1, ifelse(a$diff < -6, 2, 0))
a$diff <- ifelse(a$V4 <1,0,a$diff)
a <- a[order(a$diff),]  # 按照diff列进行升序排序
ggplot(a, aes(x = V2, y = V4, color = factor(diff), size = factor(diff))) +
  geom_point() +
  scale_color_manual(values = c("0" ="grey", "1" = "red", "2" = "grey")) +
  scale_size_manual(values = c("0" = 0.5, "1" = 0.5, "2" = 0.5))
b <- subset(a,a$diff==1)
write.csv(b, "m6Am-genes.csv")
library(ggplot2)
library(readr)
c <- runif(464,1,8540)
d <- subset(a,a$V2>3)
dim(d)
randome_matrix <- d[c,]
hist(randome_matrix$V2)
hist(b$V2)
write.csv(randome_matrix,"randome-selected-genes.csv")
plot_data <- read.table("./log2-noCap-m7G-m6Am.tsv")
plot_data <-subset(plot_data,V2 +V3 +V4 >0 )
hist(plot_data$V2)
cor(plot_data$V2,plot_data$V5)
cor(plot_data$V3,plot_data$V5)
cor(plot_data$V4,plot_data$V5)
hist(plot_data$V2 - plot_data$V3)

a<-subset(plot_data,plot_data$V2 - plot_data$V3 < 1.5)
c<-subset(plot_data,plot_data$V2 - plot_data$V3 > 3)
hist(a$V2)
hist(c$V2)
plot(a$V2,a$V5)
plot(a$V3,a$V5)
hist(a$V2)
cor(a$V2,a$V5)
cor(a$V3,a$V5)
cor(a$V4,a$V5)
plot(c$V2,c$V5)
plot(c$V3,c$V5)
cor(c$V2,c$V5)
cor(c$V3,c$V5)
hist(a$V2-a$V5)
hist(c$V2-c$V5)

b<-subset(a,a$V2 - a$V3 > 1.4 )
plot(b$V2,b$V5)
plot(b$V3,b$V5)
cor(b$V2,b$V5)
cor(b$V3,b$V5)
cor(b$V4,b$V5)
cc <- plot_data$V5-plot_data$V2
dd <- plot_data$V3-plot_data$V2
kk <- plot_data$V4-plot_data$V2
gg <- plot_data$V4-plot_data$V2
plot(cc,dd,pch=19,col="darkblue",cex=0.1)

plot(cc,kk)
cor(cc,dd)
cor(cc,kk)

plot(cc,plot_data$V4)
cor(dd,plot_data$V2)
plot(dd,plot_data$V3)
plot(cc,plot_data$V3)
cor(cc,plot_data$V3)

plot(cc,plot_data$V3)
cor(cc,plot_data$V4)

plot(plot_data$V4,plot_data$V3)


############################self-gv-bulk
a<- read.table("self-noCap-Cap-m6Acap-Am-m6Am-trans.txt")
rownames(a) <- a$V1
a$diff<-a$V3-a$V2
a$transdiff<-a$V7-a$V2
a$part <- "part"
part1<-subset(a,a$diff<=0)
barplot(part1$transdiff)
part2<-subset(a,a$diff>0)
part2<-subset(part2,part2$diff<=2)
part3<-subset(a,a$diff>2)
part3<-subset(part3,part3$diff<=4)
part4<-subset(a,a$diff>4)
a[rownames(part1),]$part <-"part1"
a[rownames(part2),]$part <-"part2"
a[rownames(part3),]$part <-"part3"
a[rownames(part4),]$part <-"part4"
a<- a[!is.na(a$part),]
cor(a$transdiff,a$diff)
m6Am_along<- read.table("./m6Am_expression.txt")
cor(m6Am_along$V1,m6Am_along$V2)
a$m6Am <- "grey"
m6Am <- read.table("./final-high-m6Am-genes.txt")
a[m6Am$V1,]$m6Am <- "red"
cor(a$V3-a$V2,a$V7-a$V2)
a <- a[order(a$m6Am),]  # 按照diff列进行升序排序
plot(a$V3-a$V2,a$V7-a$V2,pch=19,col=a$m6Am,cex=0.5)
library(ggplot2)
library(dplyr)
ggplot(a, aes(x= part, y = transdiff, fill = m6Am)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 0.2, alpha = 0.3) +
  geom_violin(trim = FALSE) +theme_minimal() + 
  theme(panel.grid = element_blank(),axis.line = element_line(linetype = "solid"))

p_values <- a %>%
  group_by(part) %>%
  summarize(p_value = wilcox.test(transdiff ~ m6Am)$p.value)
#####################################################
#suplymentary
abline(h=2,lwd=2,lty=4)
abline(h=-1,lwd=2,lty=4)
abline(v=2,lwd=2,lty=4)
abline(v=0.5,lwd=2,lty=4)
c<- subset(a,a$V3-a$V2 >2)
d<- subset(c,c$V7-c$V2 >2)
write.csv(d,"./Cap-INcrease-Trans.csv")
hist(a$diff)
hist(a$V6-a$V5)
a$m6Am <- a$V6-a$V5

Cap_m6Am_genes <- subset(a,a$diff >2)
hist(Cap_m6Am_genes$transdiff)
Cap_m6Am_genes1 <- subset(Cap_m6Am_genes,Cap_m6Am_genes$m6Am >0)
Cap_m6Am_genes1 <- subset(Cap_m6Am_genes1,Cap_m6Am_genes1$V5>2)
hist(Cap_m6Am_genes1$transdiff)



noCapgenes <- subset(a,a$V3-a$V2 <0.5)
write.csv(noCapgenes,"./Nocap-genes.csv")
cor(noCapgenes$V3-noCapgenes$V2,noCapgenes$V7-noCapgenes$V2)
noCapgenes <- subset(a,a$V3-a$V2 <1)
noCapgenes <- subset(noCapgenes,noCapgenes$V2>2)
noCapgenes <- subset(noCapgenes,noCapgenes$V3<0.5)

write.csv(noCapgenes,"Top-Noncap-genes.csv")
cor(noCapgenes$V3-noCapgenes$V2,noCapgenes$V7-noCapgenes$V2)
noCapgenes <- subset(a,a$V3-a$V2 <2)
cor(noCapgenes$V3-noCapgenes$V2,noCapgenes$V7-noCapgenes$V2)
noCapgenes <- subset(a,a$V3-a$V2 <3)
cor(noCapgenes$V3-noCapgenes$V2,noCapgenes$V7-noCapgenes$V2)
noCapgenes <- subset(a,a$V3-a$V2 <4)
cor(noCapgenes$V3-noCapgenes$V2,noCapgenes$V7-noCapgenes$V2)
noCapgenes <- subset(a,a$V3-a$V2 >=4)
cor(noCapgenes$V3-noCapgenes$V2,noCapgenes$V7-noCapgenes$V2)

aa <- cbind(noCapgenes$diff,noCapgenes$transdiff)
heatmap(aa)




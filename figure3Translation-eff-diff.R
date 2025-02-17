library(monocle3)
library(magrittr)
library()
setwd("~/Desktop/C2T-seq/Figure3-maternal-RNA/BH1-3")
plot_data_1 <- read.table("./trans-effciency-classes.tsv")
sorted_indices <- order(plot_data_1$V4)
plot_data_1 <- plot_data_1[sorted_indices, ]
plot(plot_data_1$V2,plot_data_1$V3,col = plot_data_1$V4 + 1L,cex=0.2,pch=20)
plot_data <- plot_data_1
plot_data <- as.data.frame(cbind(plot_data$V2,plot_data$V3,plot_data$V4))
rownames(plot_data) <- plot_data_1$V1
colnames(plot_data) <-c("GV_TE","MII_TE","Types","diff_TE")
plot_data$difTE <- plot_data$`MII-TE` - plot_data$`GV-TE`

plot_TanslationData <- read.table("./GV-MII-translation-log-nonCap.tsv")
plot(plot_TanslationData$V2,plot_TanslationData$V3,col = plot_TanslationData$V4 + 1L,cex=0.2,pch=20)
plot(plot_TanslationData$V2-plot_TanslationData$V5,plot_TanslationData$V3-plot_TanslationData$V6,col = plot_TanslationData$V4 + 1L,cex=0.2,pch=20)

setwd("~/Desktop/C2T-seq/Figure3-maternal-RNA/BH1-3")
log2_NonCap <- read.table("./log2-NonCap_GV_MII_IN.tsv",row.names = 1)
log2_Cap <- read.table("./log2-CPM-Cap_GV_MII_IN.tsv",row.names = 1)
log2_tail <- read.table("./log2-CPM-Tail_GV_MII_IN.tsv",row.names = 1)
row_sums <- rowSums(log2_NonCap)
log2_NonCap <- log2_NonCap[row_sums > 1, ]
a<- rownames(log2_NonCap)
log2_Cap <-log2_Cap[a,]
log2_tail <-log2_tail[a,]


library(ggplot2)
library(dplyr)
library(tidyr)

ggplot(data =plot_data , aes(x = factor(Types), y = plot_data$`MII-TE`)) + 
  geom_violin()

library(ggplot2)
library(tidyr)
library(dplyr)

# 假设plot_data是您的数据框，已经包含Types, GV_TE, 和 MII_TE这些列

# 将数据从宽格式转换为长格式
plot_data_long <- plot_data %>%
  gather(key="Condition", value="Value", GV_TE, MII_TE) %>%
  mutate(Types = as.factor(Types)) %>%
  arrange(Types, Condition)

# 创建一个新的顺序变量
plot_data_long$TypeCondition <- with(plot_data_long, factor(paste(Types, Condition, sep="_"),
                                                            levels = paste(rep(1:max(as.numeric(Types)), each=2),
                                                                           rep(c("GV_TE", "MII_TE"), max(as.numeric(Types))),
                                                                           sep="_")))

# 绘制箱线图和点图
p <- ggplot(plot_data_long, aes(x=TypeCondition, y=Value, fill=Condition)) +
  geom_boxplot(position=position_dodge(0.8)) +
 # geom_jitter(color="black", position=position_dodge(0.8), size=1) +
  labs(x='类型', y='表达量') +
  theme_minimal() +
  scale_fill_manual(values=c("GV_TE"="red", "MII_TE"="green"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.25, "cm")) 

# 打印图形
print(p)


plot(plot_data$GV_TE,plot_data$MII_TE)

umap_data <- read.csv("./GV-MII-gene-classes.csv",row.names = 1)
umap_data$difTE <- plot_data$diff_TE
umap_data$GVTE <- plot_data$GV_TE
umap_data$MIITE <- plot_data$MII_TE
colnames(umap_data) <- c("UMAP1","UMAP2","CLASS","DIFTE")
p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = MIITE)) +
  geom_point(size = 0.05) +
  scale_color_gradient2(low = "blue", high = "red", 
                        midpoint = 0, limits = c(-1, 1), oob = scales::squish) +
  theme_minimal() +
  labs(color = 'DIFTE') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))
print(p)


p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = GVTE)) +
  geom_point(size = 0.05) +
  scale_color_gradient2(low = "blue", high = "red", 
                        midpoint = 0, limits = c(-1, 1), oob = scales::squish) +
  theme_minimal() +
  labs(color = 'DIFTE') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))
print(p)

hist(plot_data$diff_TE)
TE_up <- subset(plot_data,plot_data$diff_TE >2)
TE_down <- subset(plot_data,plot_data$diff_TE < -3)
TE_down_Cap <- log2_Cap[rownames(TE_down),]
TE_down_tail <- log2_tail[rownames(TE_down),]
TE_down_nonCap <- log2_NonCap[rownames(TE_down),]

TE_up_Cap <- log2_Cap[rownames(TE_up),]
TE_up_tail <- log2_tail[rownames(TE_up),]
TE_up_nonCap <- log2_NonCap[rownames(TE_up),]

plot(TE_down_Cap$V2 - TE_down_nonCap$V2 ,TE_down_Cap$V3 - TE_down_nonCap$V3, cex =0.5,pch =20)
plot(TE_down_tail$V2 - TE_down_nonCap$V2,TE_down_tail$V3 - TE_down_nonCap$V3)

plot(TE_up_Cap$V2 - TE_up_nonCap$V2 ,TE_up_Cap$V3 - TE_up_nonCap$V3, cex = 0.5,pch =20)
plot(TE_up_tail$V2 - TE_up_nonCap$V2,TE_up_tail$V3 - TE_up_nonCap$V3)

TE_down_Cap$Cap_level_diff <- TE_down_Cap$V3 - TE_down_nonCap$V3 - TE_down_Cap$V2 + TE_down_nonCap$V2
TE_up_Cap$Cap_level_diff <- TE_up_Cap$V3 - TE_up_nonCap$V3 - TE_up_Cap$V2 + TE_up_nonCap$V2

TE_down_tail$tail_level_diff <- TE_down_tail$V3 - TE_down_nonCap$V3 - TE_down_tail$V2 + TE_down_nonCap$V2
TE_up_tail$tail_level_diff <- TE_up_tail$V3 - TE_up_nonCap$V3 - TE_up_tail$V2 + TE_up_nonCap$V2

############
#计算不通翻译趋势的polyA和Cap 的变化，修改down,up, cap,tail并记录
a<- subset(TE_down_tail,TE_down_tail$tail_level_diff > 0.5)

dim(a)

## TE-up: Cap_up 53; Cap_down 31
## TE-down: Cap_up 84; Cap_down 186

## TE-up: Tail_up 202; tail_down 141 
## TE-down: tail_up 154; tail_down 610 

#####################


translation_data <- read.table("./log2-RPKM-merge.tsv")
GV_MII_down <- subset(translation_data,translation_data$V2-translation_data$V3 >1)
GV_MII_up <- subset(translation_data,translation_data$V3-translation_data$V2 >1)
rownames(GV_MII_down) <- GV_MII_down$V1
rownames(GV_MII_up) <- GV_MII_up$V1
GV_MII_down <- GV_MII_down[rownames(GV_MII_down)%in%rownames(log2_NonCap),]
GV_MII_up <- GV_MII_up[rownames(GV_MII_up)%in%rownames(log2_NonCap),]
GV_MII_down_nonCap <- log2_NonCap[rownames(GV_MII_down),]
GV_MII_down_Cap <- log2_Cap[rownames(GV_MII_down),]
GV_MII_down_tail <- log2_tail[rownames(GV_MII_down),]
GV_MII_up_nonCap <- log2_NonCap[rownames(GV_MII_up),]
GV_MII_up_Cap <- log2_Cap[rownames(GV_MII_up),]
GV_MII_up_tail <- log2_tail[rownames(GV_MII_up),]

write.csv(plot_data,"./GV-MII-TE-classes.csv")






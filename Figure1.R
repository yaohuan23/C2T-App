# Alalysis the polyAome from total RNA-seq
library(ggplot2)
library(readr)
setwd("~/Desktop/C2T-seq/figure1/")

#log_ncRNA_OverZero_matrix <- read_delim("log-ncRNA-OverZero.matrix.csv", delim = "\t", escape_double = FALSE,trim_ws = TRUE)
log_RPKM_matrix <- read_delim("mESC-diff-methods.txt", delim = "\t", escape_double = FALSE,trim_ws = TRUE)
caculate_matrix <- subset(log_RPKM_matrix,log_RPKM_matrix$Total >2)

a<-log_RPKM_matrix
a$diff <- a$Tail-a$Total
non_polyA_genes <- subset(a,a$diff < -1 & a$Total >2)

a$diff <- ifelse(a$diff > 4.5, 1, ifelse(a$diff < -1, 2, 0))
a$diff <- ifelse(grepl("^Hist", a$geneName), 3, a$diff)
a$diff <- ifelse(grepl("^H2b", a$geneName), 3, a$diff)
a$diff <- ifelse(grepl("^H2a", a$geneName), 3, a$diff)
a$diff <- ifelse(grepl("^H3", a$geneName), 3, a$diff)
a$diff <- ifelse(grepl("^H4", a$geneName), 3, a$diff)
a <- a[order(a$diff),]  # 按照diff列进行升序排序
ggplot(a, aes(x = Total, y = Tail, color = factor(diff), size = factor(diff))) +
  geom_point() +
  scale_color_manual(values = c("0" = "black", "1" = "black", "2" = "darkgreen", "3" = "red")) +
  scale_size_manual(values = c("0" = 0.5, "1" = 0.5, "2" = 0.5, "3" = 1.5))


ggplot(a, aes(x = Total, y = Cap, color = factor(diff), size = factor(diff))) +
  geom_point() +
  scale_color_manual(values = c("0" = "black", "1" = "black", "2" = "darkgreen", "3" = "red")) +
  scale_size_manual(values = c("0" = 0.5, "1" = 0.5, "2" = 0.5, "3" = 1.5))


ggplot(a, aes(x = Total, y = polyASelected, color = factor(diff), size = factor(diff))) +
  geom_point() +
  scale_color_manual(values = c("0" = "black", "1" = "black", "2" = "darkgreen", "3" = "red")) +
  scale_size_manual(values = c("0" = 0.5, "1" = 0.5, "2" = 0.5, "3" = 1.5))
write.csv(non_polyA_genes,"non-polyA-genes.csv")


a<-log_RPKM_matrix
a$diff <- a$Cap-a$Total
a$diff <- ifelse(a$diff > 4.5, 1, ifelse(a$diff < -2, 2, 0))
a$diff <- ifelse(grepl("^Hist", a$geneName), 3, a$diff)
a$diff <- ifelse(grepl("^H2b", a$geneName), 3, a$diff)
a$diff <- ifelse(grepl("^H2a", a$geneName), 3, a$diff)
a$diff <- ifelse(grepl("^H3", a$geneName), 3, a$diff)
a$diff <- ifelse(grepl("^H4", a$geneName), 3, a$diff)
a <- a[order(a$diff),]  # 按照diff列进行升序排序
ggplot(a, aes(x = Total, y = Tail, color = factor(diff), size = factor(diff))) +
  geom_point() +
  scale_color_manual(values = c("0" = "black", "1" = "black", "2" = "darkgreen", "3" = "red")) +
  scale_size_manual(values = c("0" = 0.1, "1" = 0.5, "2" = 0.5, "3" = 1.5))


ggplot(a, aes(x = Total, y = Cap, color = factor(diff), size = factor(diff))) +
  geom_point() +
  scale_color_manual(values = c("0" = "black", "1" = "black", "2" = "darkgreen", "3" = "red")) +
  scale_size_manual(values = c("0" = 0.1, "1" = 0.5, "2" = 0.5, "3" = 1.5))


ggplot(a, aes(x = Total, y = polyASelected, color = factor(diff), size = factor(diff))) +
  geom_point() +
  scale_color_manual(values = c("0" = "black", "1" = "black", "2" = "darkgreen", "3" = "red")) +
  scale_size_manual(values = c("0" = 0.5, "1" = 0.5, "2" = 0.5, "3" = 1.5))

non_m7G_genes <- subset(a,a$diff >1 & a$Total >2)
write.csv(non_polyA_genes,"non-polyA-genes.csv")



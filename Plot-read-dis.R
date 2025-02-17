library(ggplot2)
library(reshape2)
Cap <- read.table("figure1/YH-mesc-mRNA-input_4G_TotalRelativedis.txt",header = FALSE)
PaOH <- read.table("figure1/YH-mesc-mRNA-input_3G_TotalRelativedis.txt",header = FALSE)
Tail <- read.table("figure1/tail_TotalRelativedis.txt",header = FALSE)
plot_data<- cbind(Cap,PaOH$V2,Tail$V2)
colnames(plot_data)<- c("position","Cap","P&OH","Tail")
row.names(plot_data)<-as.integer((plot_data$position-1)*100)
plot_data<-subset(plot_data[,2:4])
plot_data<-as.matrix(plot_data)
col_sums <- colSums(plot_data)
plot_data[,3]<-plot_data[,3]/15.3520
plot_data[,2]<-plot_data[,2]/1.6748
plot_data[,1]<-plot_data[,1]/1.6534

data_long <- melt(plot_data, id.vars = rownames(plot_data))
ggplot(data_long, aes(x = Var1, y = value, group = Var2 )) +
  geom_line()

write.csv(plot_data,"figure1-dis.csv")

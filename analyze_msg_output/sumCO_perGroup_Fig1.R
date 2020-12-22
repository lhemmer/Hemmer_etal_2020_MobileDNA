####################################################################################
#### Create figure 1 showing distribution of crossover (CO) number in sample groups
####################################################################################

#### load libraries

require(reshape2)
library(ggplot2)
library(ggforce)
library(ggsignif)
library(ggpubr)


#### load data

brks.sum <- read.csv("/Users/lhemmer/Documents/Github/Hemmer_etal_2020_MobileDNA/data/CO_sum_per_parent_FileS8.csv", header=TRUE,check.names=FALSE)


#### prepare data table 

## make a column called "class" to make groups for non-dysgenic, dysgenic with high fecundity, dysgenic with low fecundity

brks.sum$Class <- paste(brks.sum$dys.nondys, brks.sum$high.low.fecund, sep="_")
brks.sum$Class <- factor(brks.sum$Class, levels = c("dys_low","dys_high", "nondys_NA"))

## for second plot, looking at CO sum for each dysgenic parent, creating second data frame to add additional factors

brks.high <- brks.sum[brks.sum$dys.nondys=="dys" & brks.sum$high.low.fecund=="high",]
brks.high$parent <- factor(brks.high$parent, levels = c("701","4016", "4029","5011","5019","5027","5096"))


#### plot 

## function to quickly plot the mean and standard deviation

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


## plot dys high / low and non-dys first

p1 <- ggplot(brks.sum, aes(x=Class, y=CO.sum, fill=Class))+
  geom_violin(trim=T)+
  scale_fill_manual(values=c('red4', 'orangered1','blue3'))+
  labs(y ="Total CO per Individual")+
  labs(x ="")+
  theme_bw()+
  theme(legend.position = 'none')+
  scale_x_discrete(labels= c("Low Fecund Dys", "High Fecund Dys","Non-Dys"))+
  stat_summary(fun.data=data_summary,color="black",size=0.8)
#facet_zoom(x = class =="dys_low")


## plot high fecund dysgenic by parent

p2 <- ggplot(brks.high, aes(x=parent, y=CO.sum))+
  geom_violin(trim=T,fill="orangered1")+
  stat_summary(fun.data=data_summary, color="black",size=0.8)+
  labs(y ="Total CO per Individual",x="Parent")+
  theme_bw()+
  geom_signif(comparisons = list(c("701","4029")),map_signif_level = T)



#### output plot to pdf in arranged fasion

pdf("Figure1.pdf")
ggarrange(p1, p2, nrow=2, labels = c("A","B"))
dev.off()






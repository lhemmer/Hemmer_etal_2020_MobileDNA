####################################################################################
#### Create figure, bar plot of Crossover numbers per chromosome categorized between dysgenic high / low fecundity, non-dysgenic
####################################################################################

#### Load libraries

library(RColorBrewer)
library(ggpubr)

#### Import data

brks <- read.csv("/Users/lhemmer/Documents/Github/Hemmer_etal_2020_MobileDNA/data/All_breakpoints_FileS7.csv", header=TRUE)


#### set contig lengths

A_length <- 32431276
B_length <- 30397359
C_length <- 27291581
D_length <- 26722153
E_length <- 37623929
F_length <- 2019360


#### Remove F's for time being

brks <- brks[brks$contig != "F",]


#### get number of each class we are including in the figure

#dys.l <- length(unique(brks$sample[brks$dys == "dys"]))
nondys.l <- length(unique(brks$sample[brks$dys == "nondys"]))
dysL.l <- length(unique(brks$sample[brks$dys == "dys" & brks$fecund == "low"]))
dysN.l <- length(unique(brks$sample[brks$dys == "dys" & brks$fecund == "high"]))


## contigs in order, A = X, E = 2nd, D = 3rd, B = 4th, C = 5th
contigs <- c("A","E","D","B","C")



#### getting the number of crossovers catagorized by chromosome


## the full collection of COs for each category

#dys.all <- brks[brks$dys == "dys",]
non.all <- brks[brks$dys == "nondys",]
low.all <- brks[brks$dys == "dys" & brks$fecund == "low",]
high.all <- brks[brks$dys == "dys" & brks$fecund == "high",]


## function to get a summary of crossover numbers

co.sum.loop <- function(df, chroms) {
  
  # make output matrix
  df.out <- data.frame(matrix(, nrow=length(chroms), ncol=length(chroms)))
  
  # loop through chromosomes
  for (i in 1:length(chroms)) {
    tmp1 <- df[df$contig==chroms[i],]
    tmp2 <- split(tmp1$CO.pos, tmp1$sample)
    # for chromosomes with more than one crossover
    tmp2 <- lapply(tmp2, sort)
    tmp2 <- lapply(tmp2, as.numeric)
    # for chromosomes with one or zero crossovers
    tmp3 <- tmp2[which(sapply(tmp2, length) == 1)]
    
    # can fix up for number of crossovers if needs to be changed
    df.out[i,] <- c(length(tmp3[tmp3 == 0]), length(tmp3[tmp3 != 0]), length(tmp2[which(sapply(tmp2, length) == 2)]),
                         length(tmp2[which(sapply(tmp2, length) == 3)]), length(tmp2[which(sapply(tmp2, length) >= 4)]))
    rownames(df.out)[i] <- contigs[i]
  }
  colnames(df.out) <- c("0","1","2","3",">=4")
  rownames(df.out) <- c("X", "2", "3", "4", "5")
  return(as.matrix(df.out))
}


#### call function

non.all.num <- co.sum.loop(non.all, contigs)
low.all.num <- co.sum.loop(low.all, contigs)
high.all.num <- co.sum.loop(high.all, contigs)


#### get confidence intervals for our total crossover numbers

co.CI <- function(df, iter, chroms) {
  if (iter < 100) {
    return(print("Bootstrap number needs to be higher"))
  }
  else {
  ci <- rep(list(vector()),length(chroms))
  for (j in 1:iter) {
    all.num <- data.frame(matrix(, nrow=5, ncol=5))
    for (i in 1:length(contigs)) {
      tmp <- dys.all[dys.all$contig==contigs[i],]
      tmp2 <- split(tmp$breaks, tmp$sample)
      tmp2 <- lapply(tmp2,sort)
      tmp2 <- lapply(tmp2,as.numeric)
      tmp1 <- tmp2[which(sapply(tmp2, length) == 1)]
      samp <- sample(tmp2,553,replace=T)
      samp1 <- samp[which(sapply(samp, length) == 1)]
      all.num[i,] <- c(length(samp1[samp1 == 0]), length(samp1[samp1 != 0]), length(samp[which(sapply(samp, length) == 2)]),
                       length(samp[which(sapply(samp, length) == 3)]), length(samp[which(sapply(samp, length) >= 4)]))
    }
    sum.cols <- colSums(all.num)
    ci[[1]][j] <- sum.cols[1]
    ci[[2]][j] <- sum.cols[2]
    ci[[3]][j] <- sum.cols[3]
    ci[[4]][j] <- sum.cols[4]
    ci[[5]][j] <- sum.cols[5]
  }
  ci1 <- lapply(ci,sort)
  ci.out <- c(ci1[[1]][round(0.025*iter)],ci1[[1]][round(0.975*iter)],ci1[[2]][round(0.025*iter)],ci1[[2]][round(0.975*iter)],ci1[[3]][round(0.025*iter)],ci1[[3]][round(0.975*iter)],ci1[[4]][round(0.025*iter)],ci1[[4]][round(0.975*iter)],ci1[[5]][round(0.025*iter)],ci1[[5]][round(0.975*iter)])
  return(ci.out)
  }
}


#### confidence intervals for all three categories

ci.n <- co.CI(non.all, 1000, contigs)
#ci.n <- c(241,274,496,539,404,444,159,188,29,42)
ci.l <- co.CI(low.all, 1000, contigs)
#ci.l <- c(236,270,589,634,469,510,174,201,39,55)
ci.h <- co.CI(high.all, 1000, contigs)
#ci.h <- c(183,212,470,510,354,391,135,160,27,40)


#### #### #### #### #### #### #### #### 

#### making longform data for ggplot, it's a little messy but we want the confidence intervals near the top of the bars

#### melt data and turn it into proportional data, add confidence intervals

non.melt <- melt((non.all.num/(length(unique(non.all$sample))*5)))

high.melt <- melt((high.all.num/(length(unique(high.all$sample))*5)))

low.melt <- melt((low.all.num/(length(unique(low.all$sample))*5)))

## turn confidence intervals into proportions

ci.n1 <- ci.n/(length(unique(non.all$sample))*5)
ci.h1 <- ci.h/(length(unique(high.all$sample))*5)
ci.l1 <- ci.n/(length(unique(low.all$sample))*5)

## combine

all.sum <- rbind(high.melt,low.melt,non.melt)

## add values manually

all.sum[,4] <- c(rep("High",25),rep("Low",25),rep("Non",25))

## once again, only want CI on the very end of the bars
all.sum[,5] <- c(ci.h1[1],rep(NA,4),ci.h1[3],rep(NA,4),ci.h1[5],rep(NA,4),ci.h1[7],rep(NA,4),ci.h1[9],rep(NA,4),ci.h1[1],rep(NA,4),ci.l1[3],rep(NA,4),ci.l1[5],rep(NA,4),ci.l1[7],rep(NA,4),ci.l1[9],rep(NA,4),ci.n1[1],rep(NA,4),ci.n1[3],rep(NA,4),ci.n1[5],rep(NA,4),ci.n1[7],rep(NA,4),ci.n1[9],rep(NA,4))
all.sum[,6] <- c(ci.h1[2],rep(NA,4),ci.h1[4],rep(NA,4),ci.h1[6],rep(NA,4),ci.h1[8],rep(NA,4),ci.h1[10],rep(NA,4),ci.l1[2],rep(NA,4),ci.l1[4],rep(NA,4),ci.l1[6],rep(NA,4),ci.l1[8],rep(NA,4),ci.l1[10],rep(NA,4),ci.n1[2],rep(NA,4),ci.n1[4],rep(NA,4),ci.n1[6],rep(NA,4),ci.n1[8],rep(NA,4),ci.n1[10],rep(NA,4))

colnames(all.sum) <- c("chromosome","CO","number","dys","low","high")

#### plot

pdf("Figure2.pdf")
ggplot(data=all.sum, aes(x=dys, y=number, fill=chromosome)) + 
  geom_bar(stat="identity" ) + 
  facet_grid(~CO)+
  geom_errorbar(aes(ymin=low, ymax=high), width=0.5, position="identity")+
  labs(title="CO Count by Chromosome", y = "Proportion of Chromosomes", x = "Dysgenic or Non-Dysgenic Class") +
  scale_fill_grey(start = 0.8, end = 0.2,name="Chromosome")+
  #scale_fill_brewer(palette="Dark2") +
  #scale_fill_manual(values=c( "#810F7C" ,"#88419D", "#8C6BB1", "#8C96C6" ,"#9EBCDA")) +
  theme_minimal()+
  #theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold",hjust = 0.5))+
  scale_y_continuous(limits = c(0,0.45), expand = c(0, 0))

dev.off()







ggplot(data=all.dys, aes(x=dys, y=number, fill=chromosome)) + 
  geom_bar(stat="identity" ) + 
  facet_grid(~CO)+
  geom_errorbar(aes(ymin=low, ymax=high), width=0.5, position="identity")+
  labs(title="CO Count by Chromosome", y = "Proportion of Chromosomes", x = "Dysgenic or Non-Dysgenic Class") +
  scale_fill_grey(start = 0.8, end = 0.2,name="Chromosome")+
  #scale_fill_brewer(palette="Dark2") +
  #scale_fill_manual(values=c( "#810F7C" ,"#88419D", "#8C6BB1", "#8C96C6" ,"#9EBCDA")) +
  theme_minimal()+
  #theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold",hjust = 0.5))+
  scale_y_continuous(limits = c(0,0.45), expand = c(0, 0))





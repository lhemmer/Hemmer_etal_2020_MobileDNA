####################################################################################
#### Showing recombination rates along the length of each chromosome in Figure 3
####################################################################################


#### Import data

brks.sum <- read.csv("/Users/lhemmer/Documents/Github/Hemmer_etal_2020_MobileDNA/data/All_breakpoints_FileS7.csv", header=TRUE,check.names=FALSE)

options(scipen=999)

#### process data

## remove F contigs
brks <- brks[brks$contig != "F",]

#make parent a character rather than numeric
brks$parent <- as.character(brks$parent)

## remove all non-recombined chromosomes
brks1 <- brks[brks$breaks != 0,]


## load chromosome lengths

A_length <- 32431276
B_length <- 30397359
C_length <- 27291581
D_length <- 26722153
E_length <- 37623929
#F_length <- 2019360


#### create a function to quickly calculate recombination rate in centiMorgans per megabase

recomb.window <- function(co.vect, uniq.sample.length, chrom.length, window.size) {
  
  ## establish sliding window
  winds <- seq(from=window.size/2,to=chrom.length,by=window.size)
  window=window.size/2
  
  ## vector of recombination frequencies
  coe=rep(0,length(winds))
  
  ## we want to output in terms of centiMorgans per megabase, get adjuster
  adj.cM <- 1
  if (window.size < 1000000) {
    adj.cM <- 1000000 / window.size
  }
  
  ## calculate centimorgans per megabase within a window
  for(i in 1:length(winds)) {
    coe[i]=(adj.cM * 100) * (length(co.vect[co.vect>=winds[i]-window & co.vect<=winds[i]+window]) / uniq.sample.length)}
  
  return(coe)
}


#### just a small function for looping through categories of samples

dys.loop <- function(df, contig, contig.size, window.size) {
  
  ## separate out part of dataframe with particular contig
  df.tmp <- df[df$contig==contig,]

  ## list for output
  dys.list <- list()

  dys.list[[1]] <- df.tmp[df.tmp$fecund=="high" & df.tmp$dys=="dys",]
  dys.list[[2]] <- df.tmp[df.tmp$fecund=="low" & df.tmp$dys=="dys",]
  dys.list[[3]] <- df.tmp[df.tmp$dys=="dys",]
  dys.list[[4]] <- df.tmp[df.tmp$dys=="nondys",]
  
  ## loop recomb.window function through the list
  dys.list.out <- list()

  for (i in 1:length(dys.list)) {
    dys.list.out[[i]] <- recomb.window(dys.list[[i]]$breaks, length(unique(dys.list[[i]]$sample)), contig.size, window.size)
  }
  
  names(dys.list.out) <- c("dys.high","dys.low","dys.all","nondys")

  ## output
  return(dys.list.out)
}


#### function to consolidate and melt data for plotting in ggplot

dys.ggplot.prep <- function(df.list, contig, contig.size, window.size, cutoff) {
  
  ## run function inside of this look to limit redundancy, or just import your own dataset
  #df.list.processed <- df.list
  df.list.processed <- dys.loop(df = df.list, contig = contig, contig.size = contig.size, window.size = 500000)
  
  ## get x axis genomic position for the sliding window used in recomb.window function
  winds <- seq(from=window.size/2,to=contig.size,by=window.size)
  
  ## combine data needed for the graph, for right now only high fecund dys, low fecund dys, and non-dysgenic, add column names
  tmp.df1 <- cbind(winds,df.list.processed$dys.high,df.list.processed$dys.low,df.list.processed$nondys)
  colnames(tmp.df1) <- c("pos","high","low","non")
  tmp.df1 <- as.data.frame(tmp.df1)
  
  ## add row of empty values for the beginning of the chromosome, then melt dataframe
  tmp.df1 <- rbind(tmp.df1,c(0,0,0,0))
  tmp.df2 <- melt(tmp.df1,id=c("pos"),variable_name = "dys",value.name="rec")
  
  ## order dataframe
  tmp.df3 <- tmp.df2[order(tmp.df2$pos),]
  rownames(tmp.df3) <- c(1:length(tmp.df3$pos))
  
  ## shorten data frame so that we remove extra zero recombination regions of the genome that will shift the curve fit of recombination
  tmp.df3 <- tmp.df3[1:cutoff,]
  colnames(tmp.df3) <- c("pos","dys","rec")
  
  ## export dataframe for ggplot
  return(tmp.df3)
}


#### get ggplot data for each chromosome and store plot to variable

## X chromosome

cmX <- dys.ggplot.prep(df.list = brks1, contig = "A", contig.size = A_length, window.size = 500000, cutoff = 180)

px <- ggplot(cmX, aes(x=pos, y=rec, color=dys, shape=dys))+
  labs(title="X Chromosome",x="Genomic Position (Mb)",y="Recombination Rate (cM / Mb)")+
  geom_point(size=2) + 
  geom_smooth(method=loess, se=TRUE, fullrange=FALSE, aes(fill=dys))+
  geom_vline(xintercept=c(29250000),linetype='dashed')+ ## showing the extent of centromere suppresion of recombination
  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c('orangered1','red4', 'blue3'))+
  scale_fill_manual(values=c('orangered1','red4', 'blue3'))+
  scale_y_continuous(limits=c(-2,16))+
  scale_x_continuous(limits=c(0,A_length),breaks=seq(0,A_length,by=5000000),labels = seq(0,A_length/1000000,by=5))+
  coord_cartesian(ylim=c(0, 16))+
  theme_bw()+
  theme(legend.position="none", axis.line = element_line(size = 1), panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), panel.border = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))


## second chromosome

cm2 <- dys.ggplot.prep(df.list = brks1, contig = "E", contig.size = E_length, window.size = 500000, cutoff = 225)

p2 <- ggplot(cm2, aes(x=pos, y=rec, color=dys, shape=dys))+
  labs(title="2nd Chromosome",x="Genomic Position (Mb)",y="Recombination Rate (cM / Mb)")+
  geom_point(size=2) + 
  geom_smooth(method=loess, se=TRUE, fullrange=FALSE, aes(fill=dys))+
  geom_vline(xintercept=c(36750000),linetype='dashed')+
  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c('orangered1','red4', 'blue3'))+
  scale_fill_manual(values=c('orangered1','red4', 'blue3'))+
  scale_y_continuous(limits=c(-5,16))+
  scale_x_continuous(limits=c(0,E_length),breaks=seq(0,E_length,by=5000000),labels = seq(0,E_length/1000000,by=5))+
  coord_cartesian(ylim=c(0, 16))+
  theme_bw()+
  theme(legend.position="none", 
        axis.line = element_line(size = 1), panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), panel.border = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))


## third chromosome 

cm3 <- dys.ggplot.prep(df.list = brks1, contig = "D", contig.size = D_length, window.size = 500000, cutoff = 159)

p3 <- ggplot(cm3, aes(x=pos, y=rec, color=dys, shape=dys))+
  labs(title="3rd Chromosome",x="Genomic Position (Mb)",y="Recombination Rate (cM / Mb)")+
  geom_point(size=2) + 
  geom_smooth(method=loess, se=TRUE, fullrange=FALSE, aes(fill=dys))+
  geom_vline(xintercept=c(25750000),linetype='dashed')+
  scale_shape_manual(values=c(15, 16, 17), labels = c("High Fecund Dys", "Low Fecund Dys", "Non-Dys"))+ 
  scale_color_manual(values=c('orangered1','red4', 'blue3'), labels = c("High Fecund Dys", "Low Fecund Dys", "Non-Dys"))+
  scale_fill_manual(values=c('orangered1','red4', 'blue3'), labels = c("High Fecund Dys", "Low Fecund Dys", "Non-Dys"))+
  scale_y_continuous(limits=c(-4,30))+
  scale_x_continuous(limits=c(0,D_length),breaks=seq(0,D_length,by=5000000),labels = seq(0,D_length/1000000,by=5))+
  coord_cartesian(ylim=c(0, 30))+
  theme_bw()+
  theme(legend.position=c(0.5,0.90), legend.direction = "horizontal", legend.title = element_blank(), legend.background = element_rect(size=3),
        axis.line = element_line(size = 1), panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), panel.border = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))


## fourth chromosome

cm4 <- dys.ggplot.prep(df.list = brks1, contig = "B", contig.size = B_length, window.size = 500000, cutoff = 174)

p4 <- ggplot(cm4, aes(x=pos, y=rec, color=dys, shape=dys))+
  labs(title="4th Chromosome",x="Genomic Position (Mb)",y="Recombination Rate (cM / Mb)")+
  geom_point(size=2) + 
  geom_smooth(method=loess, se=TRUE, fullrange=FALSE, aes(fill=dys))+
  geom_vline(xintercept=c(28250000),linetype='dashed')+
  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c('orangered1','red4', 'blue3'))+
  scale_fill_manual(values=c('orangered1','red4', 'blue3'))+
  scale_y_continuous(limits=c(-4,16),expand = c(0, 0))+
  scale_x_continuous(limits=c(0,B_length),breaks=seq(0,B_length,by=5000000),labels = seq(0,B_length/1000000,by=5))+
  coord_cartesian(ylim=c(0, 16))+
  theme_bw()+
  theme(legend.position="none", 
        axis.line = element_line(size = 1), panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), panel.border = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))


## fifth chromosome

cm5 <- dys.ggplot.prep(df.list = brks1, contig = "C", contig.size = C_length, window.size = 500000, cutoff = 156)

p5 <- ggplot(cm5, aes(x=pos, y=rec, color=dys, shape=dys))+
  #ggtitle('X Chromosome')+
  labs(title="5th Chromosome",x="Genomic Position (Mb)",y="Recombination Rate (cM / Mb)")+
  geom_point(size=2) + 
  #geom_jitter() +
  geom_smooth(method=loess, se=TRUE, fullrange=FALSE, aes(fill=dys))+
  #geom_smooth(method=lm, formula=formula se=T, fullrange=FALSE, aes(fill=dys))+
  geom_vline(xintercept=c(25250000),linetype='dashed')+
  scale_shape_manual(values=c(15, 16, 17))+ 
  scale_color_manual(values=c('orangered1','red4', 'blue3'))+
  scale_fill_manual(values=c('orangered1','red4', 'blue3'))+
  #coord_cartesian(xlim=c(0,A_length), ylim=c(0,16))+
  scale_y_continuous(limits=c(-10,16))+
  scale_x_continuous(limits=c(0,C_length),breaks=seq(0,C_length,by=5000000),labels = seq(0,C_length/1000000,by=5))+
  coord_cartesian(ylim=c(0, 16))+
  theme_bw()+
  #theme_classic()+
  theme(legend.position="none", 
        axis.line = element_line(size = 1), panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), panel.border = element_blank(),
        plot.title = element_text(face="bold",hjust = 0.5))


#### load library and export figure

library(ggpubr)

pdf("Figure3.pdf")
ggarrange(ggarrange(px, p4, nrow=2, labels = c("","D")), ggarrange(p2, p5, nrow=2, labels = c("","E")), p3, ncol=3, labels=c("A","B","C"))
dev.off()







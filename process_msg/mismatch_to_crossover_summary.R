####################################################################################
#### Processing mismatch file modded with parental data to get a summary of all crossovers (COs)
####################################################################################


#### Import data

breaks <- read.csv("all-matchMismatch_processed_mod.csv", header=TRUE)


#### Set lengths of each contig

A_length <- 32431276
B_length <- 30397359
C_length <- 27291581
D_length <- 26722153
E_length <- 37623929
F_length <- 2019360


#### combine break data while retaining other information

brks1 <- breaks[,-5]
brks2 <- breaks[,-4]
colnames(brks1)[4] <- "breaks"
colnames(brks2)[4] <- "breaks"
brks <- rbind(brks1,brks2)


## order data

brks <- brks[order(brks$plate,brks$well,brks$contig,brks$breaks),]
row.names(brks) <- seq(1:nrow(brks))


#### remove all mis-annotated "crossovers" that are within 750 kb of each other

for (i in 2:nrow(brks)) {
  if ((brks$contig[i-1] == brks$contig[i]) && (abs(brks$breaks[i]-brks$breaks[i-1])) < 750000) {
    brks$breaks[i] <- 0
    brks$breaks[i-1] <- 0
  } 
}

#### remove all mis-annotated "crossovers" that are within a certain distance from the telomere or centromere
## determined by visual expection

brks$breaks[brks$contig == "A" & brks$breaks < 500000] <- NA
brks$breaks[brks$contig == "A" & brks$breaks > 29000000] <- NA

brks$breaks[brks$contig == "B" & brks$breaks < 500000] <- NA
brks$breaks[brks$contig == "B" & brks$breaks > 28000000] <- NA

brks$breaks[brks$contig == "C" & brks$breaks < 700000] <- NA
brks$breaks[brks$contig == "C" & brks$breaks > 25000000] <- NA

brks$breaks[brks$contig == "D" & brks$breaks < 700000] <- NA
brks$breaks[brks$contig == "D" & brks$breaks > 25250000] <- NA

brks$breaks[brks$contig == "E" & brks$breaks < 700000] <- NA
brks$breaks[brks$contig == "E" & brks$breaks > 36500000] <- NA

#### remove those eliminated by our criteria
brks <- brks[complete.cases(brks),]


####################################################################################
#### get a list of all crossover positons, chromosomes with no crossovers are labeled as "0"
####################################################################################

#### load library

library (plyr)

#### establish plate and contig IDs, could also extract from previous dataframe

plates <- c("Dys2","NonDys3","Dys15", "Dys16", "Dys17", "Dys21", "Dys22", "NonDys18", "NonDys19", "NonDys20")
#plates <- unique(brks$plate)
contigs <- c("A","B","C","D","E","F")
#contigs <- unique(brks$contig)


### organizes all breaks as crossovers, those with none are labled as such

## initial dataframe
breaks_total=data.frame()

## run loop

for (i in 1:length(plates)) {
  brks.tmp <- brks[brks$plate==plates[i],]
  tmp.df=data.frame()
  for (j in 1:length(contigs)) {
    tmp1 <- brks.tmp[brks.tmp$contig==contigs[j],]
    tmp2 <- split(tmp1$breaks,tmp1$well)
    tmp2 <- lapply(tmp2, sort)
    tmp2[which(sapply(tmp2, length) == 0)] <- 0
    tmp3 <- ldply(tmp2, data.frame)
    tmp3[,3] <- contigs[j]
    tmp.df <- rbind(tmp.df,tmp3)
  }
  tmp.df[,4] <- plates[i]
  breaks_total <- rbind(breaks_total,tmp.df)
}

## give column names
colnames(breaks_total) <- c("well","breaks","contig","plate")

## order crossovers by plate, sample, contig, and position
breaks_total <- breaks_total[order(breaks_total$plate,breaks_total$well,breaks_total$contig,breaks_total$breaks),]
row.names(breaks_total) <- seq(1:nrow(breaks_total))


#### output dataframe of crossover positions

write.csv(breaks_total, "CO_summary.csv")



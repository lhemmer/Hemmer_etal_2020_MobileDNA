####################################################################################
#### Organizing the mismatch
breaks <- read.csv("/Users/lhemmer/Documents/Virilis/DysNon2/All_all_mismatch1_3.csv", header=TRUE)

A_length <- 32431276
B_length <- 30397359
C_length <- 27291581
D_length <- 26722153
E_length <- 37623929
F_length <- 2019360

brks1 <- breaks[,-5]
brks2 <- breaks[,-4]
colnames(brks1)[4] <- "breaks"
colnames(brks2)[4] <- "breaks"
brks <- rbind(brks1,brks2)
brks <- brks[order(brks$plate,brks$well,brks$contig,brks$breaks),]
row.names(brks) <- seq(1:nrow(brks))

for (i in 2:nrow(brks)) {
  if ((brks$contig[i-1] == brks$contig[i]) && (abs(brks$breaks[i]-brks$breaks[i-1])) < 750000) {
    brks$breaks[i] <- 0
    brks$breaks[i-1] <- 0
  } 
}

#masking nonbreaks
#brks$breaks[brks$breaks < 100000] <- NA
#brks$breaks[brks$contig == "A" & brks$breaks > A_length - 800000] <- NA
brks$breaks[brks$contig == "A" & brks$breaks < 500000] <- NA
brks$breaks[brks$contig == "A" & brks$breaks > 29000000] <- NA
#brks$breaks[brks$contig == "B" & brks$breaks > B_length - 1000000] <- NA
brks$breaks[brks$contig == "B" & brks$breaks < 500000] <- NA
brks$breaks[brks$contig == "B" & brks$breaks > 28000000] <- NA
#brks$breaks[brks$contig == "C" & brks$breaks > C_length - 800000] <- NA
brks$breaks[brks$contig == "C" & brks$breaks < 700000] <- NA
brks$breaks[brks$contig == "C" & brks$breaks > 25000000] <- NA
#brks$breaks[brks$contig == "D" & brks$breaks > D_length - 800000] <- NA
brks$breaks[brks$contig == "D" & brks$breaks < 700000] <- NA
brks$breaks[brks$contig == "D" & brks$breaks > 25250000] <- NA
#brks$breaks[brks$contig == "E" & brks$breaks > E_length - 1000000] <- NA
brks$breaks[brks$contig == "E" & brks$breaks < 700000] <- NA
brks$breaks[brks$contig == "E" & brks$breaks > 36500000] <- NA

brks <- brks[complete.cases(brks),]

####################################################################################

library (plyr)

plates <- c("Dys15", "Dys16", "Dys17", "Dys21", "Dys22", "NonDys18", "NonDys19", "NonDys20")
contigs <- c("A","B","C","D","E","F")

### organizes all breaks
breaks_total=data.frame()
for (i in 1:length(plates)) {
  brksx <- brks[brks$plate==plates[i],]
  dysx_total=data.frame()
  for (j in 1:length(contigs)) {
    xy <- brksx[brksx$contig==contigs[j],]
    xy2 <- split(xy$breaks,xy$well)
    xy2 <- lapply(xy2, sort)
    xy2[which(sapply(xy2, length) == 0)] <- 0
    xy3 <- ldply(xy2, data.frame)
    xy3[,3] <- contigs[j]
    dysx_total <- rbind(dysx_total,xy3)
  }
  dysx_total[,4] <- plates[i]
  breaks_total <- rbind(breaks_total,dysx_total)
}

colnames(breaks_total) <- c("well","breaks","contig","plate")
breaks_total <- breaks_total[order(breaks_total$plate,breaks_total$well,breaks_total$contig,breaks_total$breaks),]
row.names(breaks_total) <- seq(1:nrow(breaks_total))

write.csv(breaks_total, "/Users/lhemmer/Documents/Virilis/DysNon2/CO_summary4.csv")




####################################################################################
#### Dys1 second
breaks <- read.csv("/Users/lhemmer/Documents/Virilis/DysNon1/All_all_mismatch1.csv", header=TRUE)

A_length <- 32431276
B_length <- 30397359
C_length <- 27291581
D_length <- 26722153
E_length <- 37623929
F_length <- 2019360

brks1 <- breaks[,c(1,2,3,4)]
brks2 <- breaks[,c(1,2,3,5)]
colnames(brks1)[4] <- "breaks"
colnames(brks2)[4] <- "breaks"
brks <- rbind(brks1,brks2)
brks <- brks[order(brks$plate,brks$indiv,brks$contig,brks$breaks),]
row.names(brks) <- seq(1:nrow(brks))

for (i in 2:nrow(brks)) {
  if ((brks$contig[i-1] == brks$contig[i]) && (abs(brks$breaks[i]-brks$breaks[i-1])) < 750000) {
    brks$breaks[i] <- 0
    brks$breaks[i-1] <- 0
  } 
}

#masking nonbreaks
#brks$breaks[brks$breaks < 100000] <- NA
#brks$breaks[brks$contig == "A" & brks$breaks > A_length - 800000] <- NA
brks$breaks[brks$contig == "A" & brks$breaks < 500000] <- NA
brks$breaks[brks$contig == "A" & brks$breaks > 29000000] <- NA
#brks$breaks[brks$contig == "B" & brks$breaks > B_length - 1000000] <- NA
brks$breaks[brks$contig == "B" & brks$breaks < 500000] <- NA
brks$breaks[brks$contig == "B" & brks$breaks > 28000000] <- NA
#brks$breaks[brks$contig == "C" & brks$breaks > C_length - 800000] <- NA
brks$breaks[brks$contig == "C" & brks$breaks < 700000] <- NA
brks$breaks[brks$contig == "C" & brks$breaks > 25000000] <- NA
#brks$breaks[brks$contig == "D" & brks$breaks > D_length - 800000] <- NA
brks$breaks[brks$contig == "D" & brks$breaks < 700000] <- NA
brks$breaks[brks$contig == "D" & brks$breaks > 25250000] <- NA
#brks$breaks[brks$contig == "E" & brks$breaks > E_length - 1000000] <- NA
brks$breaks[brks$contig == "E" & brks$breaks < 700000] <- NA
brks$breaks[brks$contig == "E" & brks$breaks > 36500000] <- NA

brks <- brks[complete.cases(brks),]

####################################################################################

library (plyr)

plates <- c("Dys2","NonDys3")
contigs <- c("A","B","C","D","E","F")

### organizes all breaks
breaks_total=data.frame()
for (i in 1:length(plates)) {
  brksx <- brks[brks$plate==plates[i],]
  dysx_total=data.frame()
  for (j in 1:length(contigs)) {
    xy <- brksx[brksx$contig==contigs[j],]
    xy2 <- split(xy$breaks,xy$indiv)
    xy2 <- lapply(xy2, sort)
    xy2[which(sapply(xy2, length) == 0)] <- 0
    xy3 <- ldply(xy2, data.frame)
    xy3[,3] <- contigs[j]
    dysx_total <- rbind(dysx_total,xy3)
  }
  dysx_total[,4] <- plates[i]
  breaks_total <- rbind(breaks_total,dysx_total)
}

colnames(breaks_total) <- c("indiv","breaks","contig","plate")
breaks_total <- breaks_total[order(breaks_total$plate,breaks_total$indiv,breaks_total$contig,breaks_total$breaks),]
row.names(breaks_total) <- seq(1:nrow(breaks_total))

write.csv(breaks_total, "/Users/lhemmer/Documents/Virilis/DysNon1/CO_summary1.csv")











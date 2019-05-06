library(tidyverse)
library(reshape2)
library(pracma)



# limiting PCA to GOA climate variables!
dat <- read.csv("clim vars for ordination.csv", row.names=1)

# put years into row names
rownames(dat) <- dat[,1]
dat <- dat[,-1]

# dropping AL for now as I think it's redundant with SLP!
dat <- dat[,-1]

head(dat)

# remove some extra  time series
dat <- dat[,-c(6,7,9,10,11,12,15,18)]

# check
head(dat)

# and restrict to local variables
dat <- dat[,-c(1:3)]

# check
head(dat)


# check pairwise correlations for detrended data

cor(detrend(as.matrix(dat[rownames(dat) %in% 1950:1988,])), use="p"); cor(detrend(as.matrix(dat[rownames(dat) %in% 1989:2012,])), use="p")

dat1 <- dat[rownames(dat) %in% 1950:1988,]
dat2 <- dat[rownames(dat) %in% 1989:2012,]
# and fit pca by era

pca.1 <- svd(cov(scale(dat1), use="p"))
pca.2 <- svd(cov(scale(dat2), use="p"))

# randomize svd after Jackson (first era)
x <- scale(dat1)

eigenvalue <- matrix(NA, nrow = 10000, ncol = ncol(x))
set.seed(1299)
for(i in 1:10000){
  new <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for(j in 1:ncol(x)){
    new[,j] <- sample(x[,j], size = nrow(x), replace = T)
  }
  
  temp <- svd(cov(new, use="pairwise.complete.obs"))
  eigenvalue[i,] <- temp$d
  #eigenvector1[i,] <- temp$u[,1]
  #eigenvector2[i,] <- temp$u[,2]
}

evs <- matrix(NA, nrow = ncol(x), ncol = 3)
colnames(evs) <- c("mean", "LI", "UI")
for (i in 1:ncol(x)){
  eigenvalue <- eigenvalue[order(eigenvalue[,i]),]
  evs[i,1] <- mean(eigenvalue[,i])
  evs[i,2] <- eigenvalue[250,i]
  evs[i,3] <- eigenvalue[9750,i]
}

comp.pca1 <- as.data.frame(cbind(pca.1$d, evs)) 
colnames(comp.pca1)[1] <- "eigenvalue"
comp.pca1$axis <- 1:nrow(comp.pca1)

######
# and now the second era
x <- scale(dat2)

eigenvalue <- matrix(NA, nrow = 10000, ncol = ncol(x))
#eigenvector1 <- eigenvector2 <- matrix(NA, nrow = 1000, ncol = 12)
set.seed(1298)
for(i in 1:10000){
  new <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for(j in 1:ncol(x)){
    new[,j] <- sample(x[,j], size = nrow(x), replace = T)
  }
  
  temp <- svd(cov(new, use="pairwise.complete.obs"))
  eigenvalue[i,] <- temp$d
  #eigenvector1[i,] <- temp$u[,1]
  #eigenvector2[i,] <- temp$u[,2]
}

evs <- matrix(NA, nrow = ncol(x), ncol = 3)
colnames(evs) <- c("mean", "LI", "UI")
for (i in 1:ncol(x)){
  eigenvalue <- eigenvalue[order(eigenvalue[,i]),]
  evs[i,1] <- mean(eigenvalue[,i])
  evs[i,2] <- eigenvalue[250,i]
  evs[i,3] <- eigenvalue[9750,i]
}

comp.pca2 <- as.data.frame(cbind(pca.2$d, evs)) 
colnames(comp.pca2)[1] <- "eigenvalue" 


comp.pca2$axis <- 1:nrow(comp.pca2)

# and combine!
comp.pca1$era <- "1950-1988"
comp.pca2$era <- "1989-2012"

comp.pca <- rbind(comp.pca1, comp.pca2)

ev <- cbind(pca.1$d^2/sum(pca.1$d^2), pca.2$d^2/sum(pca.2$d^2))
colnames(ev) <- c("1950-1988", "1989-2010")
rownames(ev) <- 1:nrow(ev)

ev <- melt(ev)
names(ev) <- c("axis", "period", "variance.explained")

# compare variance explained for each era
ggplot(data=ev, aes(fill=period, y=variance.explained, x=axis)) + geom_bar(position="dodge", stat="identity")

# axis1 explains 82.9% in 1950-1988, 60.3% in 1989-2012
comp.pca$axis <- as.factor(comp.pca$axis)
comp.pca$axis.num <- as.numeric(comp.pca$axis)
# determine which are interpretable
pdf("scree plot GOA climate by era NOT DETRENDED.pdf", 3,4)
ggplot(comp.pca, aes(x=axis, y=eigenvalue)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=LI, ymax=UI), colour="red", width=.5) +
  geom_line(aes(x=axis.num, y=mean), colour="red") + facet_wrap(~era) + ylab("Eigenvalue") + xlab("Principal component") + theme_gray()
dev.off()

# and save for combined figure!
clim.plot <- ggplot(comp.pca, aes(x=axis, y=eigenvalue)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=LI, ymax=UI), colour="red", width=.5) +
  geom_line(aes(x=axis.num, y=mean), colour="red") + facet_wrap(~era) + ylab("Eigenvalue") + xlab("Principal component") + theme_gray() +
  ggtitle("Environmental variability")

# now the community pca
# load population time series
dat <- read.csv("salmon and non-salmon biology mar 28.csv", row.names = 1)

# drop GOAPOP
keep <- colnames(dat)!="GOAPOP"
dat <- dat[,keep] 

# and fit pca by era
dat1 = dat[rownames(dat) %in% 1965:1988,]

dat2 = dat[rownames(dat) %in% 1989:2012,]


pca.1 <- svd(cov(scale(dat1), use="p"))
pca.2 <- svd(cov(scale(dat2), use="p"))

#bootstrap svd after Jackson (first era)
x <- scale(dat1)

eigenvalue <- matrix(NA, nrow = 10000, ncol = ncol(x))
set.seed(1299)
for(i in 1:10000){
  new <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for(j in 1:ncol(x)){
    new[,j] <- sample(x[,j], size = nrow(x), replace = T)
  }
  
  temp <- svd(cov(new, use="pairwise.complete.obs"))
  eigenvalue[i,] <- temp$d
  #eigenvector1[i,] <- temp$u[,1]
  #eigenvector2[i,] <- temp$u[,2]
}

evs <- matrix(NA, nrow = ncol(x), ncol = 3)
colnames(evs) <- c("mean", "LI", "UI")
for (i in 1:ncol(x)){
  eigenvalue <- eigenvalue[order(eigenvalue[,i]),]
  evs[i,1] <- mean(eigenvalue[,i])
  evs[i,2] <- eigenvalue[250,i]
  evs[i,3] <- eigenvalue[9750,i]
}

comp.pca1 <- as.data.frame(cbind(pca.1$d, evs)) 
colnames(comp.pca1)[1] <- "eigenvalue"
comp.pca1$axis <- 1:nrow(comp.pca1)

######
# and now the second era
x <- scale(dat2)

eigenvalue <- matrix(NA, nrow = 10000, ncol = ncol(x))
#eigenvector1 <- eigenvector2 <- matrix(NA, nrow = 1000, ncol = 12)
set.seed(1298)
for(i in 1:10000){
  new <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for(j in 1:ncol(x)){
    new[,j] <- sample(x[,j], size = nrow(x), replace = T)
  }
  
  temp <- svd(cov(new, use="pairwise.complete.obs"))
  eigenvalue[i,] <- temp$d
  #eigenvector1[i,] <- temp$u[,1]
  #eigenvector2[i,] <- temp$u[,2]
}

evs <- matrix(NA, nrow = ncol(x), ncol = 3)
colnames(evs) <- c("mean", "LI", "UI")
for (i in 1:ncol(x)){
  eigenvalue <- eigenvalue[order(eigenvalue[,i]),]
  evs[i,1] <- mean(eigenvalue[,i])
  evs[i,2] <- eigenvalue[250,i]
  evs[i,3] <- eigenvalue[9750,i]
}

comp.pca2 <- as.data.frame(cbind(pca.2$d, evs)) 
colnames(comp.pca2)[1] <- "eigenvalue" 

# comp.pca2

comp.pca2$axis <- 1:nrow(comp.pca2)

# and combine!
comp.pca1$era <- "1965-1988"
comp.pca2$era <- "1989-2012"

comp.pca <- rbind(comp.pca1, comp.pca2)

ev <- cbind(pca.1$d^2/sum(pca.1$d^2), pca.2$d^2/sum(pca.2$d^2))
colnames(ev) <- c("1965-1988", "1989-2010")
rownames(ev) <- 1:nrow(ev)

ev <- melt(ev)
names(ev) <- c("axis", "period", "variance.explained")

# compare variance explained for each era
ggplot(data=ev, aes(fill=period, y=variance.explained, x=axis)) + geom_bar(position="dodge", stat="identity")
# axis1 explains 48.2% in 1950-1988, 31.9% in 1989-2012
comp.pca$axis <- as.factor(comp.pca$axis)
comp.pca$axis.num <- as.numeric(comp.pca$axis)
# determine which are interpretable
# only plot the first 10 axes
pdf("scree plot non-detrended community variability by era.pdf", 3,4)
ggplot(filter(comp.pca, axis.num <= 10), aes(x=axis, y=eigenvalue)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=LI, ymax=UI), colour="red", width=.5) +
  geom_line(aes(x=axis.num, y=mean), colour="red") + facet_wrap(~era) + ylab("Eigenvalue") + xlab("Principal component") + theme_grey()
dev.off()
# only axis 1 in early period, none in second period!!!

# save for a combined plot

comm.plot <- ggplot(filter(comp.pca, axis.num <= 10), aes(x=axis, y=eigenvalue)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=LI, ymax=UI), colour="red", width=.5) +
  geom_line(aes(x=axis.num, y=mean), colour="red") + facet_wrap(~era) + ylab("Eigenvalue") + xlab("Principal component") + theme_grey() + 
  ggtitle("Community variability")

# combine
pdf("scree plot PCA environmental and community not detrended.pdf", 8, 5) 
ggarrange(clim.plot,  comm.plot, labels = c("a)", "b)"),  widths=c(1,1), ncol=2)
dev.off()
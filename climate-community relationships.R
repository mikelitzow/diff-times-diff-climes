library(ggplot2)
library(reshape2)
library(dplyr)
library(pracma)
library(zoo)
library(MuMIn)
library(ncdf4)
library(chron)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(broom)
library(lemon)

# begin with a plot of Aleutian Low SLPa distributions in the two eras

# load the SLP data as described in methods/SI
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/climate data/prmsl.mon.mean.7.28.15.nc")

raw <- ncvar_get(nc, "time")
h <- raw/24
d <- dates(h, origin = c(1,1,1800))

# Pick start and end dates (Jan 1950-Dec 2012):
d <- d[949:1704]  

# just the box in spatial change in SD plot! 48-56 deg. N, 192-206 deg. E:
sd.x <- ncvar_get(nc, "lon", start=97, count=8)
sd.y <- ncvar_get(nc, "lat", start=18, count=5)

SLP1 <- ncvar_get(nc, "prmsl", start=c(97,18,949), count=c(8,5,length(d)))

SLP1 <- aperm(SLP1, 3:1)  # reverse order of dimensions
SLP1 <- SLP1[,5:1,]  # Reverse order of latitudes to be increasing for convenience (in later plotting)
sd.y <- rev(sd.y)  # Also reverse corresponding vector of latitudes

SLP1 <- matrix(SLP1, nrow=dim(SLP1)[1], ncol=prod(dim(SLP1)[2:3]))  # Change to matrix
# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(sd.y, length(sd.x))   # Vector of latitudes
lon <- rep(sd.x, each = length(sd.y))   # Vector of longitudes

dimnames(SLP1) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

m <- months(d)
f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series

p.mu1 <- apply(SLP1, 2, f)	# Compute monthly means for each time series (location)
p.mu1 <- p.mu1[rep(1:12, round(length(d)/12)),] 

# Compute matrix of anomalies!
SLP1.anom <- (SLP1 - p.mu1)

# and smoothing with 11-mo rolling mean
SLP1.anom <- rollmean(SLP1.anom, 11, fill=NA)
mean.anom <- rowMeans(SLP1.anom)

# make a data frame
slp <- data.frame(slp.11 = mean.anom, year=rep(1950:2012, each=12), month=1:12)
slp$era <- ifelse(slp$year <= 1988, "1950-1988", "1989-2012")

slp$dec.yr <- slp$year+(slp$month-0.5)/12

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

AL.p <- ggplot(slp, aes(slp.11, fill=era)) + 
  geom_density(alpha=0.8) +
  xlim(-700,700) + 
  scale_fill_manual(values=cb[c(6,8)]) +
  theme_bw() +
  theme(legend.position = c(0.2,0.8), legend.title = element_blank()) +
  xlab("Sea level pressure anomaly (Pa)") +
  ggtitle("Aleutian Low SLPa distribution")

# and compare the variance between the two eras
# using a randomization test
F.rand <- NA

for(i in 1:10000){
 # i <- 1
  temp.x <- sample(slp$slp.11, length(slp$slp.11[slp$year<=1988]), replace=T)
  temp.y <- sample(slp$slp.11, length(slp$slp.11[slp$year>1988]), replace=T)
  
  F.rand[i] <- var.test(temp.x, temp.y)$statistic
}

sum(F.rand > var.test(slp$slp.11[slp$year<=1988], slp$slp.11[slp$year>1988])$statistic)

var.test(slp$slp.11[slp$year<=1988], slp$slp.11[slp$year>1988])

# F test to compare two variances
# 
# data:  slp$slp.11[slp$year <= 1988] and slp$slp.11[slp$year > 1988]
# F = 1.843, num df = 462, denom df = 282, p-value = 3.144e-08
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   1.490123 2.267081
# sample estimates:
#   ratio of variances 
# 1.842961

# load biology and climate data
b.dat <- read.csv("GOA community data.csv", row.names = 1)
c.dat <- read.csv("GOA environmental data.csv", row.names = 1)

# first, scale biology
b.dat <- scale(b.dat[rownames(b.dat) %in% 1965:2012,])

# and put into early and late climate DFs
clim.early <- as.data.frame(cbind(c.dat[rownames(c.dat) %in% 1950:1988,c(1,2,7)]))
clim.late <- as.data.frame(cbind(c.dat[rownames(c.dat) %in% 1989:2012,c(1,2,7)]))

sm.early <- rollmean(clim.early,2, align="right", fill=NA)
rownames(sm.early) <- 1950:1988

sm.late <- as.data.frame(rollmean(clim.late,2, align="right", fill=NA))
rownames(sm.late) <- 1989:2012

# limit early climate data years to match biology
clim.early <- clim.early[rownames(clim.early) >= 1965,]
sm.early <- as.data.frame(sm.early[rownames(sm.early) >= 1965,])

# now make output objects and find era-specific relationships!
pdo.e <- sst.e <- npgo.e <- pdo.l <- sst.l <- npgo.l <- NA

# also keeping track of whether raw or smoothed data were superior in each case
best.mod <- matrix(nrow=14, ncol=3)
dimnames(best.mod) <- list(colnames(b.dat), c("PDO", "SST", "NPGO"))

# now loop through each biology ts!

for(i in 1:ncol(b.dat)){
 # i <- 1
  # first, loop through the early climate-biology relationships
  pdo1 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ clim.early$PDO, na.action = "na.exclude")
  pdo2 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ sm.early$PDO, na.action = "na.exclude")
  
  ifelse(AICc(pdo1) < AICc(pdo2), pdo.e[i] <- summary(pdo1)$coefficients[2,1], 
         pdo.e[i] <- summary(pdo2)$coefficients[2,1])
  
  # record best model for later plotting
  ifelse(AICc(pdo1) < AICc(pdo2), best.mod[i,1] <- 1, best.mod[i,1] <- 2)

  sst1 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ clim.early$SST, na.action = "na.exclude")
  sst2 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ sm.early$SST, na.action = "na.exclude")
  
  ifelse(AICc(sst1) < AICc(sst2), sst.e[i] <- summary(sst1)$coefficients[2,1], 
         sst.e[i] <- summary(sst2)$coefficients[2,1])
  
  # record best model for later plotting
  ifelse(AICc(sst1) < AICc(sst2), best.mod[i,2] <- 1, best.mod[i,2] <- 2)
  
  npgo.1 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ clim.early$NPGO, na.action = "na.exclude")
  npgo.2 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ sm.early$NPGO, na.action = "na.exclude")
  
  ifelse(AICc(npgo.1) < AICc(npgo.2), npgo.e[i] <- summary(npgo.1)$coefficients[2,1], 
         npgo.e[i] <- summary(npgo.2)$coefficients[2,1])

  # record best model for later plotting
  ifelse(AICc(npgo.1) < AICc(npgo.2), best.mod[i,3] <- 1, best.mod[i,3] <- 2)
  
  # and now the late era
  # using either smoothed or unsmoothed data depending on best model for early era!
  ifelse(AICc(pdo1) < AICc(pdo2), x <- clim.late$PDO, x <- sm.late$PDO)
  mod <- lm(b.dat[rownames(b.dat) %in% 1989:2012,i] ~ x, na.action = "na.exclude")
  pdo.l[i] <- summary(mod)$coefficients[2,1]
  
  ifelse(AICc(sst1) < AICc(sst2), x <- clim.late$SST, x <- sm.late$SST)
  mod <- lm(b.dat[rownames(b.dat) %in% 1989:2012,i] ~ x, na.action = "na.exclude")
  sst.l[i] <- summary(mod)$coefficients[2,1]
  
  ifelse(AICc(npgo.1) < AICc(npgo.2), x <- clim.late$NPGO, x <- sm.late$NPGO)
  mod <- lm(b.dat[rownames(b.dat) %in% 1989:2012,i] ~ x, na.action = "na.exclude")
  npgo.l[i] <- summary(mod)$coefficients[2,1]
  }

hist(pdo.e); hist(pdo.l)
hist(sst.e); hist(sst.l)
hist(npgo.e); hist(npgo.l)

t.test(abs(pdo.e), abs(pdo.l), paired = T)
# Paired t-test
# 
# data:  abs(pdo.e) and abs(pdo.l)
# t = 3.0875, df = 13, p-value = 0.008653
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.09178231 0.51952216
# sample estimates:
#   mean of the differences 
# 0.3056522 

t.test(abs(sst.e), abs(sst.l), paired = T)
# Paired t-test

# data:  abs(sst.e) and abs(sst.l)
# t = 4.4293, df = 13, p-value = 0.00068
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.4184594 1.2153419
# sample estimates:
#   mean of the differences 
# 0.8169006 

t.test(abs(npgo.e), abs(npgo.l), paired = T)
# Paired t-test
# 
# data:  abs(npgo.e) and abs(npgo.l)
# t = 3.6067, df = 13, p-value = 0.003191
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.08482157 0.33821883
# sample estimates:
#   mean of the differences 
# 0.2115202

hist(pdo.e); hist(pdo.l)
hist(sst.e); hist(sst.l)

sc.pdo <- scale(c(pdo.e, pdo.l), center=F)[1:28]
sc.sst <- scale(c(sst.e, sst.l), center=F)[1:28]
sc.npgo <- scale(c(npgo.e, npgo.l), center=F)[1:28]

plot.cors <- data.frame(Predictor=c(rep("PDO",28), rep("SST",28), rep("NPGO",28)), 
                        Era=rep(c(rep("1965-1988",14), rep("1989-2012",14)),3), 
                        value=abs(c(pdo.e, pdo.l, sst.e, sst.l, npgo.e, npgo.l)), 
                        scaled.value=abs(c(sc.pdo, sc.sst, sc.npgo)))

plot.cors$Predictor <- reorder(plot.cors$Predictor, rep(c(1,3,2), each=28))

clim.p <- ggplot(plot.cors, aes(scaled.value)) + geom_density(aes(fill=Era), color="black", alpha=0.8) +
  facet_wrap(~ Predictor) + 
  theme_bw() + ggtitle("Community response to climate") +
  scale_fill_manual(values=cb[c(6,8)]) +
  theme(legend.position = c(0.18,0.88), legend.title = element_blank(), legend.key.size = unit(4, "mm"), legend.margin = margin(1,1,1,1, "mm")) +
  xlab("Scaled regression coefficient (absolute value)") +
  xlim(0,3) 


# calculate pairwise correlations in each era
e.cor <- cor(b.dat[rownames(b.dat) <= 1988,], use="p")

l.cor <- cor(b.dat[rownames(b.dat) > 1989,], use="p")

keep <- lower.tri(e.cor)
e.cor <- e.cor[keep]
l.cor <- l.cor[keep]

summary(abs(e.cor)); summary(abs(l.cor))

t.test(abs(e.cor), abs(l.cor), paired = T)
# Paired t-test
# 
# data:  abs(e.cor) and abs(l.cor)
# t = 6.3116, df = 90, p-value = 1.019e-08
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1315436 0.2523950
# sample estimates:
#   mean of the differences 
# 0.1919693 


# and combine for a plot
comm.cors <- data.frame(era=rep(c("1965-1988", "1989-2012"), each=length(e.cor)), cor=c(abs(e.cor), abs(l.cor)))

comm.p <- ggplot(comm.cors, aes(cor)) + geom_density(aes(fill=era), color="black", alpha=0.8) +
  xlab("Pairwise correlation (absolute value)") +
  theme_bw() + ggtitle("Pairwise community correlations") +
  scale_fill_manual(values=cb[c(6,8)]) +
  theme(legend.position = c(0.6,0.8), legend.title = element_blank()) +
  xlab("Pairwise correlation (absolute value)") +
  xlim(0,1)

pdf("AL and climate-community and community-community plot.pdf",4,10)
ggarrange(AL.p, clim.p,  comm.p, labels = c("a)", "b)", "c)"), ncol=1, nrow=3)
dev.off()


# now make SI plots of the regressions in early and late eras...

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make new list of label names
colnames(b.dat)
good.names <- c("Sablefish", "Walleye pollock", "Arrowtooth flounder", "Pacific herring", "Tanner crab",
                "Sidestripe shrimp", "Capelin", "Northern shrimp", "Pacific cod", "Chum salmon", 
                "Pink salmon", "Coho salmon", "Sockeye salmon", "Northern rockfish")

# first, PDO
pdo.plot <- data.frame(TS=NA, year=NA, Era=NA, x=NA, y=NA)

for(i in 1:ncol(b.dat)){
 # i <- 1
  temp <- data.frame(TS=colnames(b.dat)[i], year=1965:2012, 
                     Era=c(rep("1965-1988", length(1965:1988)), rep("1989-2012", length(1989:2012))), x=NA, y=b.dat[,i])
  ifelse(best.mod[i,1] == 1, temp$x <- clim.early$FMA.PDO, temp$x <- sm.early$FMA.PDO)
  
  pdo.plot <- rbind(pdo.plot, temp)
}

# drop leading NA
pdo.plot <- pdo.plot[2:nrow(pdo.plot),]

labels <- paste(good.names, " (", best.mod[,1], " yr.)", sep="")
pdo.plot$labels <- rep(labels, each=length(1965:2012))

pp <- ggplot(pdo.plot, aes(x=x, y=y, color=Era)) + geom_point() + facet_wrap(~labels) + 
  scale_color_manual(values=c("1965-1988" = cb[6], "1989-2012" = cb[8])) + 
  geom_smooth(method = "lm", se=F) + theme_bw() +
  xlab("PDO") + ylab("Standard anomaly") +
  theme(strip.text=element_text(size=9),legend.title = element_blank())

png("community responses PDO SI Fig.png", 7,6, units="in", res=300)
reposition_legend(pp, 'left', panel = 'panel-4-3')
dev.off()

pdf("community responses PDO SI Fig.pdf", 6,6)
ggplot(pdo.plot, aes(x=x, y=y, color=Era)) + geom_point() + facet_wrap(~labels) + 
  scale_fill_manual(values=c("1965-1988" = cb[6], "1989-2012" = cb[8])) + 
  geom_smooth(method = "lm", se=F) + theme_gray() + xlab("PDO") + ylab("Standard anomaly")
dev.off()

# and NPGO
npgo.plot <- data.frame(TS=NA, year=NA, Era=NA, x=NA, y=NA)

for(i in 1:ncol(b.dat)){
  # i <- 1
  temp <- data.frame(TS=colnames(b.dat)[i], year=1965:2012, 
                     Era=c(rep("1965-1988", length(1965:1988)), rep("1989-2012", length(1989:2012))), x=NA, y=b.dat[,i])
  ifelse(best.mod[i,1] == 1, temp$x <- clim.early$FMA.NPGO, temp$x <- sm.early$FMA.NPGO)
  
  npgo.plot <- rbind(npgo.plot, temp)
}

# drop leading NA
npgo.plot <- npgo.plot[2:nrow(npgo.plot),]

labels <- paste(good.names, " (", best.mod[,1], " yr.)", sep="")
npgo.plot$labels <- rep(labels, each=length(1965:2012))

pp <- ggplot(npgo.plot, aes(x=x, y=y, color=Era)) + geom_point() + facet_wrap(~labels) + 
  scale_color_manual(values=c("1965-1988" = cb[6], "1989-2012" = cb[8])) + 
  geom_smooth(method = "lm", se=F) + theme_bw() +
  xlab("NPGO") + ylab("Standard anomaly") +
  theme(strip.text=element_text(size=9),legend.title = element_blank())

png("community responses npgo SI Fig.png", 7,6, units="in", res=300)
reposition_legend(pp, 'left', panel = 'panel-4-3')
dev.off()

pdf("community responses npgo SI Fig.pdf", 6,6)
ggplot(npgo.plot, aes(x=x, y=y, color=Era)) + geom_point() + facet_wrap(~labels) + 
  scale_fill_manual(values=c("1965-1988" = cb[6], "1989-2012" = cb[8])) + 
  geom_smooth(method = "lm", se=F) + theme_gray() + xlab("npgo") + ylab("Standard anomaly")
dev.off()



# now SST
sst.plot <- data.frame(TS=NA, year=NA, Era=NA, x=NA, y=NA)

for(i in 1:ncol(b.dat)){
 #  i <- 1
  temp <- data.frame(TS=colnames(b.dat)[i], year=1965:2012, 
                     Era=c(rep("1965-1988", length(1965:1988)), rep("1989-2012", length(1989:2012))), x=NA, y=b.dat[,i])
  ifelse(best.mod[i,2] == 1, temp$x <- clim.early$FMA.SST, temp$x <- sm.early$FMA.SST)
  
  sst.plot <- rbind(sst.plot, temp)
}

# drop leading NA
sst.plot <- sst.plot[2:nrow(sst.plot),]

labels <- paste(good.names, " (", best.mod[,1], " yr.)", sep="")
sst.plot$labels <- rep(labels, each=length(1965:2012))

pdf("community responses SST SI Fig.pdf", 6,6)
ggplot(sst.plot, aes(x=x, y=y, color=Era)) + geom_point() + facet_wrap(~labels) + 
  geom_smooth(method = "lm", se=F) + theme(legend.position = "top") + xlab("SST (ºC)") + ylab("Scaled anomaly")
dev.off()


tp <- ggplot(sst.plot, aes(x=x, y=y, color=Era)) + geom_point() + facet_wrap(~labels) + 
  scale_color_manual(values=c("1965-1988" = cb[6], "1989-2012" = cb[8])) + 
  geom_smooth(method = "lm", se=F) + theme_bw() +
  xlab("SST (ºC)") + ylab("Standard anomaly") +
  theme(strip.text=element_text(size=9),legend.title = element_blank())

png("community responses SST SI Fig.png", 7,6, units="in", res=300)
reposition_legend(tp, 'left', panel = 'panel-4-3')
dev.off()

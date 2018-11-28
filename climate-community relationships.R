library(ggplot2)
library(reshape2)
library(dplyr)
library(pracma)
library(zoo)
library(MuMIn)

# load biology data
b.dat <- read.csv("salmon and non-salmon biology mar 28.csv", row.names = 1)

# drop GOAPOP
keep <- colnames(b.dat)!="GOAPOP"
b.dat <- b.dat[,keep] 

# load climate data
c.dat <- read.csv("clim vars for ordination.csv", row.names=1)

# put years into row names
rownames(c.dat) <- c.dat[,1]
c.dat <- c.dat[,-1]

# dropping AL for now as I think it's redundant with SLP!
c.dat <- c.dat[,-1]

# remove some extra  time series
c.dat <- c.dat[,-c(6,7,9,10,11,12,15,18)]

# check
head(c.dat)

# first, scale biology
b.dat <- scale(b.dat[rownames(b.dat) %in% 1965:2012,])

# now calculate era-specific PC1 of local climate variables
# NOT detrending as we are interested in effects of actual variability, not detrended variability!

dat1 <- as.matrix(c.dat[rownames(c.dat) %in% 1950:1988,4:10])
dat2 <- as.matrix(c.dat[rownames(c.dat) %in% 1989:2012,4:10])

# and fit pca by era

pca.1 <- svd(cov(scale(dat1), use="p"))
pca.2 <- svd(cov(scale(dat2), use="p"))

# remove NAs and calculate PC scores from svd
x.na <- is.na(dat1)
x.0 <- dat1
x.0[x.na==T] <- 0
pc1.early <- x.0 %*% pca.1$u[,1]

x.na <- is.na(dat2)
x.0 <- dat2
x.0[x.na==T] <- 0
pc1.late <- x.0 %*% pca.2$u[,1]

# and put into early and late climate DFs
clim.early <- as.data.frame(cbind(c.dat[rownames(c.dat) %in% 1950:1988,c(2,8)], 
                    pc1.early))
clim.late <- as.data.frame(cbind(c.dat[rownames(c.dat) %in% 1989:2012,c(2,8)], pc1.late))

sm.early <- rollmean(clim.early,2, align="right", fill=NA)
rownames(sm.early) <- 1950:1988

sm.late <- as.data.frame(rollmean(clim.late,2, align="right", fill=NA))
rownames(sm.late) <- 1989:2012

# limit early climate data years to match biology
clim.early <- clim.early[rownames(clim.early) >= 1965,]
sm.early <- as.data.frame(sm.early[rownames(sm.early) >= 1965,])

# now make output objects and find era-specific relationships!
pdo.e <- sst.e <- pc1.e <- pdo.l <- sst.l <- pc1.l <- NA

# also keeping track of whether raw or smoothed data were superior in each case
best.mod <- matrix(nrow=14, ncol=3)
dimnames(best.mod) <- list(colnames(b.dat), c("PDO", "SST", "PC1"))

# now loop through each biology ts!

for(i in 1:ncol(b.dat)){
  #i <- 1
  # first, loop through the early climate-biology relationships
  pdo1 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ clim.early$FMA.PDO, na.action = "na.exclude")
  pdo2 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ sm.early$FMA.PDO, na.action = "na.exclude")
  
  ifelse(AICc(pdo1) < AICc(pdo2), pdo.e[i] <- summary(pdo1)$coefficients[2,1], 
         pdo.e[i] <- summary(pdo2)$coefficients[2,1])
  
  # record best model for later plotting
  ifelse(AICc(pdo1) < AICc(pdo2), best.mod[i,1] <- 1, best.mod[i,1] <- 2)

  sst1 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ clim.early$FMA.SST, na.action = "na.exclude")
  sst2 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ sm.early$FMA.SST, na.action = "na.exclude")
  
  ifelse(AICc(sst1) < AICc(sst2), sst.e[i] <- summary(sst1)$coefficients[2,1], 
         sst.e[i] <- summary(sst2)$coefficients[2,1])
  
  # record best model for later plotting
  ifelse(AICc(sst1) < AICc(sst2), best.mod[i,2] <- 1, best.mod[i,2] <- 2)
  
  pc1.1 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ clim.early$pc1.early, na.action = "na.exclude")
  pc1.2 <- lm(b.dat[rownames(b.dat) %in% 1965:1988,i] ~ sm.early$pc1.early, na.action = "na.exclude")
  
  ifelse(AICc(pc1.1) < AICc(pc1.2), pc1.e[i] <- summary(pc1.1)$coefficients[2,1], 
         pc1.e[i] <- summary(pc1.2)$coefficients[2,1])

  # record best model for later plotting
  ifelse(AICc(pc1.1) < AICc(pc1.2), best.mod[i,3] <- 1, best.mod[i,3] <- 2)
  
  # and now the late era
  # using either smoothed or unsmoothed data depending on best model for early era!
  ifelse(AICc(pdo1) < AICc(pdo2), x <- clim.late$FMA.PDO, x <- sm.late$FMA.PDO)
  mod <- lm(b.dat[rownames(b.dat) %in% 1989:2012,i] ~ x, na.action = "na.exclude")
  pdo.l[i] <- summary(mod)$coefficients[2,1]
  
  ifelse(AICc(sst1) < AICc(sst2), x <- clim.late$FMA.SST, x <- sm.late$FMA.SST)
  mod <- lm(b.dat[rownames(b.dat) %in% 1989:2012,i] ~ x, na.action = "na.exclude")
  sst.l[i] <- summary(mod)$coefficients[2,1]
  
  ifelse(AICc(pc1.1) < AICc(pc1.2), x <- clim.late$pc1.late, x <- sm.late$pc1.late)
  mod <- lm(b.dat[rownames(b.dat) %in% 1989:2012,i] ~ x, na.action = "na.exclude")
  pc1.l[i] <- summary(mod)$coefficients[2,1]
  }

hist(pdo.e); hist(pdo.l)
hist(sst.e); hist(sst.l)
hist(-pc1.e); hist(-pc1.l)

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

t.test(abs(pc1.e), abs(pc1.l), paired = T)
# Paired t-test
# 
# data:  abs(pc1.e) and abs(pc1.l)
# t = 0.88656, df = 13, p-value = 0.3914
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.04605126  0.11015432
# sample estimates:
#   mean of the differences 
# 0.03205153 
hist(pdo.e); hist(pdo.l)
hist(sst.e); hist(sst.l)
hist(-pc1.e); hist(-pc1.l)

plot.cors <- data.frame(Predictor=c(rep("PDO",28), rep("SST",28), rep("Climate PC1",28)), 
                        Era=rep(c(rep("1965-1988",14), rep("1989-2012",14)),3), 
                        value=abs(c(pdo.e, pdo.l, sst.e, sst.l, -pc1.e, -pc1.l)))
# adjust bins
xb <- 7

pdf("histograms community response to climate.pdf", 6, 4)
ggplot(plot.cors, aes(value, fill=Era)) + geom_histogram(bins=xb, position = "dodge") +
  facet_wrap(~Predictor) + xlab("Regression coefficient (absolute value)")#+ geom_freqpoly(aes(color=Era), bins=xb)
dev.off()

png("histograms community response to climate.png", 6, 4, units="in", res=300)
ggplot(plot.cors, aes(value)) + geom_histogram(bins=xb, position = "dodge", fill="gray65", color="black") +
  facet_grid(Era ~ Predictor, scales="free_y") + xlab("Regression coefficient (absolute value)") +
  theme(axis.title = element_text(size=16)) + theme(axis.text = element_text(size=16)) + 
  theme(strip.text = element_text(size=16))
  
dev.off()

restr.cors <- filter(plot.cors, Predictor != "Climate PC1")
png("histograms community response to PDO and SST.png", 4, 4, units="in", res=300)
ggplot(restr.cors, aes(value)) + geom_histogram(bins=xb, position = "dodge", fill="gray65", color="black") +
  facet_grid(Era ~ Predictor, scales="free_x") + xlab("Regression coefficient (absolute value)") +
  theme(axis.title = element_text(size=14)) + theme(axis.text = element_text(size=12)) + 
  theme(strip.text = element_text(size=14))

dev.off()

# save for combined plot
xb <-7
clim.p <- ggplot(restr.cors, aes(value)) + geom_histogram(bins=xb, position = "dodge", fill="gray65", color="black") +
  facet_grid(Era ~ Predictor, scales="free_x") + xlab("Regression coefficient (absolute value)") +
  theme_gray() + ggtitle("Community response to climate")

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
xb <- 13
comm.p <- ggplot(comm.cors, aes(cor)) + geom_histogram(bins=xb, position = "dodge", fill="gray65", color="black") +
  facet_grid( ~ era) + xlab("Pairwise correlation (absolute value)") +
  theme_gray() + ggtitle("Pairwise community correlations")
comm.p

png("climate-community and community-community plot.png", 7,4, units="in", res=300)
ggarrange(clim.p,  comm.p, labels = c("a)", "b)"), ncol=2)
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
  scale_color_manual(values=c("1965-1988" = cb[6], "1989-2012" = cb[7])) + 
  geom_smooth(method = "lm", se=F) + theme_gray() +
  xlab("PDO") + ylab("Standard anomaly") +
  theme(strip.text=element_text(size=9),legend.title = element_blank())

png("community responses PDO SI Fig.png", 7,6, units="in", res=300)
reposition_legend(pp, 'left', panel = 'panel-4-3')
dev.off()

pdf("community responses PDO SI Fig.pdf", 6,6)
ggplot(pdo.plot, aes(x=x, y=y, color=Era)) + geom_point() + facet_wrap(~labels) + 
  scale_fill_manual(values=c("1965-1988" = cb[6], "1989-2012" = cb[7])) + 
  geom_smooth(method = "lm", se=F) + theme_gray() + xlab("PDO") + ylab("Standard anomaly")
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
  scale_color_manual(values=c("1965-1988" = cb[6], "1989-2012" = cb[7])) + 
  geom_smooth(method = "lm", se=F) + theme_gray() +
  xlab("SST (ºC)") + ylab("Standard anomaly") +
  theme(strip.text=element_text(size=9),legend.title = element_blank())

png("community responses SST SI Fig.png", 7,6, units="in", res=300)
reposition_legend(tp, 'left', panel = 'panel-4-3')
dev.off()
# note that the following was mistaken bc I scaled climate values by era, so the mean = 0 for each era by definition!
# # and plot distributions for each predictor in the two eras!
# names(clim.early) <- names(clim.late) <- c("PDO", "SST", "Climate PC1")
# c.e <- melt(clim.early)
# c.l <- melt(clim.late)
# 
# c.e$Era <- "1965-1988"
# c.l$Era <- "1989-2012"
# 
# plot.dist <- rbind(c.e, c.l)
# 
# # adjust bins
# xb <- 10
# png("histograms climate distribution.png", 6, 4, units="in", res=300)
# ggplot(plot.dist, aes(value)) + geom_histogram(bins=xb, position = "dodge", fill="gray65", color="black") +
#   facet_grid(Era ~ variable, scales="free") + xlab("Scaled anomaly")+
#   theme(axis.title = element_text(size=16)) + theme(axis.text = element_text(size=16)) + 
#   theme(strip.text = element_text(size=16))
# dev.off()

# # compare climate variables between eras using GLS
# plot.dist$Era <- as.factor(plot.dist$Era)
# mod <- gls(value ~ Era, data=plot.dist[plot.dist$variable == "Climate PC1",], correlation = corAR1())
# summary(mod)

# recheck with unscaled data!
c.dat$Era <- ifelse(rownames(c.dat) <= 1988, "early", "late")
mod <- gls(FMA.PDO ~ Era, data=c.dat, correlation = corAR1())
summary(mod)
# Generalized least squares fit by REML
# Model: FMA.PDO ~ Era 
# Data: c.dat 
# AIC      BIC    logLik
# 172.5112 180.9547 -82.25559
# 
# Correlation Structure: AR(1)
# Formula: ~1 
# Parameter estimate(s):
#   Phi 
# 0.5619826 
# 
# Coefficients:
#   Value Std.Error     t-value p-value
# (Intercept) -0.08671959 0.3138777 -0.27628463  0.7833
# Eralate     -0.03506383 0.4784999 -0.07327865  0.9418
# 
# Correlation: 
#   (Intr)
# Eralate -0.588
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -2.20732882 -0.57093044 -0.07342859  0.73640997  1.91964036 
# 
# Residual standard error: 1.088773 
# Degrees of freedom: 63 total; 61 residual

mod <- gls(FMA.SST ~ Era, data=c.dat, correlation = corAR1())
summary(mod)

# Generalized least squares fit by REML
# Model: FMA.SST ~ Era 
# Data: c.dat 
# AIC      BIC    logLik
# 104.8313 113.2748 -48.41567
# 
# Correlation Structure: AR(1)
# Formula: ~1 
# Parameter estimate(s):
#   Phi 
# 0.3160111 
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)  4.814915 0.1180728 40.77923  0.0000
# Eralate     -0.002896 0.1882872 -0.01538  0.9878
# 
# Correlation: 
#   (Intr)
# Eralate -0.61 
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -2.40144423 -0.59530800  0.06699666  0.74432154  2.05446626 
# 
# Residual standard error: 0.5391645 
# Degrees of freedom: 63 total; 61 residual



#######
# now climate-biology residuals!
regr.dat <- as.data.frame(b.dat)
# add era
regr.dat$era <- "early "
regr.dat$era[rownames(regr.dat) >= 1989] <- "late"

regr.dat$PDO <- c.dat$FMA.PDO[match(rownames(regr.dat), rownames(c.dat))]
regr.dat$SST <- c.dat$FMA.SST[match(rownames(regr.dat), rownames(c.dat))]

# now calculate PC1 for local climate over entire time series!
pca <- svd(cov(scale(c.dat[,4:10]), use="p"))

# remove NAs and calculate PC scores from svd
x.na <- is.na(c.dat[,4:10])
x.0 <- scale(c.dat[,4:10])
x.0[x.na==T] <- 0
pc1 <- x.0 %*% pca$u[,1]

regr.dat$clim.PC1 <- pc1[match(rownames(regr.dat), rownames(pc1))]

# and smooth each climate variable

regr.dat$PDO.sm <- rollmean(regr.dat$PDO, 2, align="right", fill=NA)
regr.dat$SST.sm <- rollmean(regr.dat$SST, 2, align="right", fill=NA)
regr.dat$clim.PC1.sm <- rollmean(regr.dat$clim.PC1, 2, align="right", fill=NA)


# and objects for residuals in stationary and non-stationary models
pdo.st <- sst.st <- pc1.st <- pdo.nst <- sst.nst <- pc1.nst <- as.data.frame(matrix(ncol=2, nrow=ncol(b.dat)))
colnames(pdo.st) <- colnames(sst.st) <- colnames(pc1.st) <- colnames(pdo.nst) <- colnames(sst.nst) <- 
  colnames(pc1.nst) <- c("AR(1)", "P")
# now loop through each biology ts!
library(car)

for(i in 1:ncol(b.dat)){
  #i <- 1
  # first, loop through the early climate-biology relationships
  pdo1 <- lm(regr.dat[rownames(regr.dat) %in% 1965:1988,i] ~ regr.dat$PDO[rownames(regr.dat) %in% 1965:1988], na.action = "na.exclude")
  pdo2 <- lm(regr.dat[rownames(regr.dat) %in% 1965:1988,i] ~ regr.dat$PDO.sm[rownames(regr.dat) %in% 1965:1988], na.action = "na.exclude")
  
  ifelse(AICc(pdo1) < AICc(pdo2), x <- regr.dat$PDO,  x <- regr.dat$PDO.sm)
  
  pdo.st[i,1] <- dwt(lm(regr.dat[,i] ~ x, na.action="na.omit"))$r
  pdo.st[i,2] <- dwt(lm(regr.dat[,i] ~ x, na.action="na.omit"))$p
  
  # now non-stationary model!
  pdo.nst[i,1] <- dwt(lm(regr.dat[,i] ~ x*regr.dat$era, na.action="na.omit"))$r
  pdo.nst[i,2] <- dwt(lm(regr.dat[,i] ~ x*regr.dat$era, na.action="na.omit"))$p
  
  # same for sst
  sst1 <- lm(regr.dat[rownames(regr.dat) %in% 1965:1988,i] ~ regr.dat$SST[rownames(regr.dat) %in% 1965:1988], na.action = "na.exclude")
  sst2 <- lm(regr.dat[rownames(regr.dat) %in% 1965:1988,i] ~ regr.dat$SST.sm[rownames(regr.dat) %in% 1965:1988], na.action = "na.exclude")
  
  ifelse(AICc(sst1) < AICc(sst2), x <- regr.dat$SST,  x <- regr.dat$SST.sm)
  
  sst.st[i,1] <- dwt(lm(regr.dat[,i] ~ x, na.action="na.omit"))$r
  sst.st[i,2] <- dwt(lm(regr.dat[,i] ~ x, na.action="na.omit"))$p
  
  # now non-stationary model!
  sst.nst[i,1] <- dwt(lm(regr.dat[,i] ~ x*regr.dat$era, na.action="na.omit"))$r
  sst.nst[i,2] <- dwt(lm(regr.dat[,i] ~ x*regr.dat$era, na.action="na.omit"))$p 
  
  # and for climate PC1
  pc1.1 <- lm(regr.dat[rownames(regr.dat) %in% 1965:1988,i] ~ regr.dat$clim.PC1[rownames(regr.dat) %in% 1965:1988], na.action = "na.exclude")
  pc1.2 <- lm(regr.dat[rownames(regr.dat) %in% 1965:1988,i] ~ regr.dat$clim.PC1.sm[rownames(regr.dat) %in% 1965:1988], na.action = "na.exclude")
  
  ifelse(AICc(pc1.1) < AICc(pc1.2), x <- regr.dat$clim.PC1,  x <- regr.dat$clim.PC1.sm)
  
  pc1.st[i,1] <- dwt(lm(regr.dat[,i] ~ x, na.action="na.omit"))$r
  pc1.st[i,2] <- dwt(lm(regr.dat[,i] ~ x, na.action="na.omit"))$p
  
  # now non-stationary model!
  pc1.nst[i,1] <- dwt(lm(regr.dat[,i] ~ x*regr.dat$era, na.action="na.omit"))$r
  pc1.nst[i,2] <- dwt(lm(regr.dat[,i] ~ x*regr.dat$era, na.action="na.omit"))$p
}

# combine and plot

pdo.st$type <- sst.st$type <- pc1.st$type <- "Stationary"
pdo.nst$type <- sst.nst$type <- pc1.nst$type <- "Non-stationary"

pdo.st$Predictor <- pdo.nst$Predictor <- "PDO"
sst.st$Predictor <- sst.nst$Predictor <- "SST"
pc1.st$Predictor <- pc1.nst$Predictor <- "Climate PC1"

resid.plot <- rbind(pdo.st, pdo.nst, sst.st, sst.nst, pc1.st, pc1.nst)
names(resid.plot)[1] <- "AR"
#resid.plot <- melt(resid.plot, measure.vars = c("AR(1)", "P"))

# adjust bins
xb <- 10
pdf("histograms climate-biology residuals", 6, 4)
ggplot(resid.plot, aes(AR)) + geom_histogram(bins=xb, position = "dodge") +
  facet_grid(type ~ Predictor, scales="free") #+ xlab("Scaled anomaly")#+ geom_freqpoly(aes(color=Era), bins=xb)
dev.off()
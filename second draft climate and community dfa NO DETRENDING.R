library(ggplot2)
library(reshape2)
library(dplyr)
library(pracma)
library(MARSS)
library(tidyr)
# devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(cowplot)
library(gridExtra)
library(broom)
library(lemon)
library(MuMIn)

# load data...
dat <- read.csv("clim vars for ordination.csv", row.names=1)

# put years into row names
rownames(dat) <- dat[,1]
dat <- dat[,-1]

# dropping AL and NPI as I think we will leave out these dynamics to keep the paper simple
# the relevant atmosphere-ocean dynamics have been addequately addressed in PRSB ms. for our purposes!
dat <- dat[,-c(1,2)]

head(dat)

# remove some extra  time series
dat <- dat[,-c(5,6,8,9,10,11,14,17)]


# change names to plot-friendly labels
colnames(dat) <- c("PDO", "NPGO", "SLP gradient", "Freshwater", "Wind stress", "Downwelling", "SST", "Advection", "SSH")

# plot TS for SI
do.dat <- as.data.frame(scale(dat)) # scale to plot on 1 axis
plot.dat <- gather(do.dat)
plot.dat$year <- 1950:2012

# make an order key for plotting
plot.dat$order <- rep(c(9,8,4,3,7,2,5,1,4), each=nrow(dat))

# and reorder labels
plot.dat$key <- reorder(plot.dat$key, plot.dat$order)
# save as SI Fig

pdf("SI - climate time series.pdf", 8,6)
ggplot(plot.dat, aes(x=year, y=value)) + geom_bar(position="dodge", stat="identity") + facet_wrap(~key) +
  ylab("Standard anomaly") + xlab("") + theme_gray() + geom_hline(yintercept = 0, col="dark grey")
dev.off()


# and remove PDO/NPGO
dat <- dat[,-c(1,2)]

# check
head(dat)
# now fit DFA models
# first the early era

# limit to PDO era
e.cli.dat = as.matrix(dat[rownames(dat) %in% 1950:1988,])

# and transpose
e.cli.dat <- t(e.cli.dat)

# refit the different error structures and compare

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()
# fit lots of models & store results
# NOTE: this will take a long time to run!
for(R in levels.R) {
  for(m in 1:3) {  # allowing up to 3 trends
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(e.cli.dat, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data

# 2 trend equal var covar the best!
write.csv(model.data, "climate dfa table 1950-1988 NOT DETRENDED data.csv")
early.model.data <- model.data

############
# and the late!
# limit to NPGO era 
l.cli.dat = as.matrix(dat[rownames(dat) %in% 1989:2012,])

# and transpose
l.cli.dat <- t(l.cli.dat)
model.data = data.frame()
# fit lots of models & store results
# NOTE: this will take a long time to run!
for(R in levels.R) {
  for(m in 1:3) {
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(l.dat, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data
write.csv(model.data, "climate dfa table 1989-2012 NOT DETRENDED data.csv")

##
# fit best model for each era

model.list = list(A="zero", m=2, R="equalvarcov") # best model for early era
mod.early = MARSS(e.cli.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# need to rotate!
# and rotate the loadings
Z.est = coef(mod.early, type="matrix")$Z
H.inv = varimax(coef(mod.early, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H.inv)

# and plot
# reverse trend 2 to plot
Z.rot[,2] <- -Z.rot[,2]

Z.rot$names <- rownames(e.cli.dat)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(e.cli.dat)
Z.rot$plot.names <- reorder(Z.rot$names, 1:14)
ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar(stat="identity", position="dodge") +
  theme_gray() 
# the two trends are largely mirror images of each other. Trend 1 relates SST to upwelling, advection,
# FW discharge, and the SLP gradient

# and the late era

# equal var covar one trend is the best for this era!
model.list = list(A="zero", m=1, R="equalvarcov")
mod.late = MARSS(l.cli.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...
modCI.late <- MARSSparamCIs(mod.late)
modCI.late

# aside...get correlations among loadings!
cor.test(modCI.late$par$Z, Z.rot[,1])
# Pearson's product-moment correlation
# 
# data:  modCI.late$par$Z and Z.rot[, 1]
# t = -0.8767, df = 5, p-value = 0.4208
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.8770059  0.5351367
# sample estimates:
# cor 
# -0.3650208

cor.test(modCI.late$par$Z, Z.rot[,2])
# Pearson's product-moment correlation
# 
# data:  modCI.late$par$Z and Z.rot[, 2]
# t = 0.95548, df = 5, p-value = 0.3832
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.5114694  0.8843198
# sample estimates:
# cor 
# 0.3929334 
plot.CI.late <- data.frame(names=rownames(l.cli.dat), mean=modCI.late$par$Z, upCI=modCI.late$par.upCI$Z,
                           lowCI=modCI.late$par.lowCI$Z)

plot.CI.late <- arrange(plot.CI.late, mean)
plot.CI.late$names.order <- reorder(plot.CI.late$names, plot.CI.late$mean)
dodge <- position_dodge(width=0.9)

ggplot(plot.CI.late, aes(x=names.order, y=mean)) + geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme(axis.text.x = element_text(angle = 90, vjust=0, hjust=1)) +
  geom_hline(yintercept = 0) + theme_gray()

# best model returns noise!

# combine into one plot
Z.rot$key <- rep(c("Trend 1", "Trend 2"), each=7)
plot.CI.late$names.order <- reorder(plot.CI.late$names, plot.CI.late$mean)

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

early.plot <- ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar(stat="identity", position="dodge") +
  theme_gray() + ylab("Loading") + xlab("") + ggtitle("1950-1988") + 
  scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3])) +
  theme(legend.position = c(0.8,0.2), legend.title=element_blank()) + geom_hline(yintercept = 0) +
  theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=12)) 

# reorder Z plot
plot.CI.late$plot.names <- reorder(plot.CI.late$names, c(1,3,5,2,6,7,4))

late.plot <- ggplot(plot.CI.late, aes(x=plot.names, y=mean)) + geom_bar(position="dodge", stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme_gray() + theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=12)) +
  geom_hline(yintercept = 0) + ggtitle("1989-2012")

# combine
pdf("loadings plot detrended climate dfa both eras.pdf", 6, 4)
ggarrange(early.plot,  late.plot, labels = c("a)", "b)"),  widths=c(1.3,1), ncol=2)
dev.off()


###
# now compare with best 1-trend model for the early era
# equal var covar one trend is the best for this era!
model.list = list(A="zero", m=1, R="equalvarcov")
mod.early.1 = MARSS(e.cli.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...
modCI.early.1<- MARSSparamCIs(mod.early.1)
modCI.early.1

plot.CI.early.1 <- data.frame(names=rownames(e.cli.dat), mean=modCI.early.1$par$Z, upCI=modCI.early.1$par.upCI$Z,
                           lowCI=modCI.early.1$par.lowCI$Z)

# plot.CI.early.1 <- arrange(plot.CI.early.1, mean)
plot.CI.early.1$names.order <- reorder(plot.CI.late$names, plot.CI.late$mean)
dodge <- position_dodge(width=0.9)

ggplot(plot.CI.early.1, aes(x=names.order, y=mean)) + geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme(axis.text.x = element_text(angle = 90, vjust=0, hjust=1)) +
  geom_hline(yintercept = 0) + theme_gray()

# I'm not sure that makes sense as a model!
# Paired t-test
#
# data:  plot.CI.early$mean and plot.CI.late$mean
# t = 3.4201, df = 6, p-value = 0.01414
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.04861737 0.29310266
# sample estimates:
#   mean of the differences
# 0.17086

# plot loadings together
plot.loads <- rbind(plot.CI.early, plot.CI.late)
plot.loads$era <- c(rep("1965-1988", 7), rep("1989-2012", 7))

pdf("one trend loadings climate NOT DETRENDED.pdf", 4,3.5)
ggplot(plot.loads, aes(x=names.order, y=mean, fill=era)) + theme_grey() + geom_bar(position="dodge", stat="identity") +
  geom_hline(yintercept = 0, color="dark grey")+ geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) + ylab("Loading") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust=0, hjust=1), legend.position = c(0.18, 0.85), legend.title = element_blank())

dev.off()


# now the whole time series

dat = as.matrix(dat)

# and transpose
cli.dat <- t(dat)


# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()
# # fit lots of models & store results
# # NOTE: this will take a long time to run!
for(R in levels.R) {
  for(m in 1:3) {
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(cli.dat, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

model.data$d.AICc <- model.data$AICc-min(model.data$AICc)
model.data <- arrange(model.data, d.AICc)
model.data

write.csv(model.data, "dfa table climate NOT detrended all years.csv")

###########
# now...the community analysis
dat <- read.csv("salmon and non-salmon biology mar 28.csv", row.names = 1)

# drop GOAPOP
keep <- colnames(dat)!="GOAPOP"
dat <- dat[,keep] 

# check
head(dat)

# limit to 1965:2012
dat <- dat[rownames(dat) <= 2012,]

# change names to plot-friendly labels
colnames(dat) <- c("Sablefish", "Walleye pollock", "Arrowtooth flounder", "Pacific herring", "Tanner crab", "Sidestripe shrimp", "Capelin", "Northern shrimp", "Pacific cod", "Chum salmon", "Pink salmon",
                   "Coho salmon", "Sockeye salmon", "Northern rockfish")

# plot TS for SI
do.dat <- as.data.frame(scale(dat)) # scale to plot on 1 axis
plot.dat <- gather(do.dat)
plot.dat$year <- 1965:2012

# save as SI Fig

pdf("SI - community time series.pdf", 7,6)
ggplot(plot.dat, aes(x=year, y=value)) + geom_bar(position="dodge", stat="identity") + facet_wrap(~key) +
  ylab("Standard anomaly") + xlab("") + geom_hline(yintercept = 0, col="dark grey") + theme_gray()
dev.off()

# fit to the early era
e.comm.dat = as.matrix(dat[rownames(dat) %in% 1965:1988,])

# and transpose
e.comm.dat <- t(e.comm.dat)
# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()
# fit lots of models & store results
# NOTE: this will take a long time to run!
for(R in levels.R) {
  for(m in 1:3) {
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(e.comm.dat, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

arrange(model.data, AICc)
#save
write.csv(arrange(model.data, AICc), "DFA model selection NOT DETRENDED community 65-88 IMPROVED CONVERGENCE.csv")

early.model.data <- model.data

# and fit to the late era

l.comm.dat = as.matrix(dat[rownames(dat) %in% 1989:2012,])

# and transpose
l.comm.dat <- t(l.comm.dat)

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()
# fit lots of models & store results
# NOTE: this will take a long time to run!
for(R in levels.R) {
  for(m in 1:3) {
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(l.comm.dat, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

arrange(model.data, AICc)
#save
write.csv(arrange(model.data, AICc), "DFA model selection NOT DETRENDED community 89-12 IMPROVED CONVERGENCE.csv")

# #####
# # compare the loadings for each best community model
# #####
##
# fit best model for each era

model.list = list(A="zero", m=2, R="diagonal and unequal") # best model for early era
mod.early = MARSS(e.comm.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# need to rotate!
# and rotate the loadings
Z.est = coef(mod.early, type="matrix")$Z
H.inv = varimax(coef(mod.early, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H.inv)

# and plot
# reverse trend 2 to plot
#Z.rot[,2] <- -Z.rot[,2]

Z.rot$names <- rownames(e.comm.dat)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(e.comm.dat)
Z.rot$plot.names <- reorder(Z.rot$names, 1:28)
ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar(stat="identity", position="dodge") +
  theme_gray() 

# and rotate the trends!
trends.rot = solve(H.inv) %*% mod.early$states

# plot these trends!
plot.trends <- as.data.frame(t(trends.rot)) %>%
  gather()
plot.trends$year <- 1965:1988

ggplot(plot.trends, aes(year, value, color=key)) + geom_line() + theme_gray()

# and the late era

# equal var covar one trend is the best for this era!
model.list = list(A="zero", m=1, R="diagonal and equal")
mod.late = MARSS(l.comm.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...
modCI.late <- MARSSparamCIs(mod.late)
modCI.late

plot.CI.late <- data.frame(names=rownames(l.comm.dat), mean=modCI.late$par$Z, upCI=modCI.late$par.upCI$Z,
                           lowCI=modCI.late$par.lowCI$Z)

plot.CI.late <- arrange(plot.CI.late, mean)
plot.CI.late$names.order <- reorder(plot.CI.late$names, c(13,12,9,7,2,1,11,10,14,3,4,5,8,6))
dodge <- position_dodge(width=0.9)

ggplot(plot.CI.late, aes(x=names.order, y=mean)) + geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme(axis.text.x = element_text(angle = 90, vjust=0, hjust=1)) +
  geom_hline(yintercept = 0) + theme_gray()


# combine into one plot
Z.rot$key <- rep(c("Trend 1", "Trend 2"), each=14)
#plot.CI.late$names.order <- reorder(plot.CI.late$names, plot.CI.late$mean)

early.plot <- ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar( position="dodge",stat="identity") +
  theme_gray() + ylab("Loading") + xlab("") + ggtitle("1965-1988") + 
  scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3])) +
  theme(legend.position = c(0.12,0.85), legend.title=element_blank()) + geom_hline(yintercept = 0) +
  theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=12)) 

late.plot <- ggplot(plot.CI.late, aes(x=names.order, y=mean)) + 
  geom_bar(position="dodge", stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme_gray() + theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=12)) +
  geom_hline(yintercept = 0) + ggtitle("1989-2012")

# combine
pdf("loadings plot non-detrended community dfa both eras.pdf", 8, 5)
ggarrange(early.plot,  late.plot, labels = c("a)", "b)"),  widths=c(1.3,1), ncol=2)
dev.off()

#############################################
# and test for significant AR(1) patterns
# now calculate difference between predicted (trend * loading) and observed for each trend!
# and rotate the loadings
Z.est = coef(mod.early, type="matrix")$Z
H.inv = varimax(coef(mod.early, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H.inv)

Z.rot$names <- rownames(e.comm.dat)

# and rotate the trends!
trends.rot = solve(H.inv) %*% mod.early$states

t1.r <- t2.r  <- matrix(nrow=24, ncol=14)
# t.r = matrix(nrow=48, ncol=14)
# pred = matrix(nrow=48, ncol=14)
for(i in 1:14){
  # i <- 1
  t1.r[,i] <- (trends.rot[1,] * Z.rot[i,1]) - e.comm.dat[i,]
  t2.r[,i] <- (trends.rot[2,] * Z.rot[i,2]) - e.comm.dat[i,]

  #  t4.r[,i] <- (trends.rot[4,] * Z.rot[i,4]) - all.dat[i,]
  # pred[,i] <- (as.matrix(Z.rot[,1:4]) %*% trends.rot)[i,]
  # t.r[,i] <- pred[,i] - all.dat[i,]
}

colnames(t1.r) <- colnames(t2.r) <- colnames(t3.r) <- colnames(t4.r) <- rownames(all.comm.dat)

# colnames(t.r) <- rownames(all.dat)


# calculate AR(1) values for the resulting time series
# using Durbin-Watson test
library(lmtest)
dw.1 <- dw.2 <- NA

for(i in 1:14){
  dw.1[i] <- dwtest(t1.r[,i] ~ 1)$p.value
  dw.2[i] <- dwtest(t2.r[,i] ~ 1)$p.value
  # dw.4[i] <- dwtest(t4.r[,i] ~ 1)$p.value
}

sum(dw.1 <=0.05)
sum(dw.2 <=0.05)

# and fit AR(1) models to each
ar.1 <- ar.2 <- NA

for(i in 1:14){
  ar.1[i] <- ar(na.omit(t1.r[,i]), order.max=1, aic=F)$ar
  ar.2[i] <- ar(na.omit(t2.r[,i]), order.max=1, aic=F)$ar
  # dw.4[i] <- dwtest(t4.r[,i] ~ 1)$p.value
}

# and save for later comparison!
comm.65.88.ar <- data.frame(ar=c(ar.1, ar.2))

# now ar for residuals in 2nd era
late.res <- residuals(mod.late)$residuals[1:14,]

# replace 0s with NAs
drop <- late.res==0
late.res[drop] <- NA

ar.1 <- NA
for(i in 1:14){
  ar.1[i] <- ar(na.omit(t1.r[,i]), order.max=1, aic=F)$ar
}

# plot, out of curiosity!
late.res <- gather(as.data.frame(t(late.res)))
late.res$year <- 1989:2012
ggplot(late.res, aes(year, value)) + geom_line() + geom_point() + facet_wrap(~key) + theme_gray()

# and save
comm.89.12.ar <- data.frame(ar=ar.1)

##############################################
# and calculate the correlations among loadings
cor.test(modCI.late$par$Z, Z.rot[,1])
# Pearson's product-moment correlation
# 
# data:  modCI.late$par$Z and Z.rot[, 1]
# t = -0.092438, df = 12, p-value = 0.9279
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.5494776  0.5111388
# sample estimates:
# cor 
# -0.02667495 

cor.test(modCI.late$par$Z, Z.rot[,2])
# Pearson's product-moment correlation
# 
# data:  modCI.late$par$Z and Z.rot[, 2]
# t = -2.1364, df = 12, p-value = 0.05393
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.82557475  0.00782567
# sample estimates:
# cor 
# -0.5249334 

###############################################
# now find the best model for entire community time series
all.comm.dat = dat

# and transpose
all.comm.dat <- t(all.comm.dat)

cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()
res.ar <- matrix(nrow=1, ncol=17)
# fit lots of models & store results
# NOTE: this will take a long time to run!
for(R in levels.R) {
  for(m in 1:3) {
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(all.comm.dat, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
    
    # add in calculation of residual AR(1) values!
     dw.p <- NA
     res <- residuals(kemz)$residuals
    for(ii in 1:nrow(res)){
      
      dw.p[ii] <- dwtest(res[ii,] ~ 1)$p.value
     
    }
    
res.ar <- rbind(res.ar, dw.p)
    
  } # end m loop
} # end R loop

res.ar
colnames(res.ar) <- rownames(res)
res.ar <- na.omit(res.ar)

model.res.table <- as.data.frame(res.ar)
model.res.table$R <- model.data$R
model.res.table$m <- model.data$m

write.csv(model.res.table, "residual ar1 all comm models 1965-2012.csv")

#res.ar$R <- model.data$R
###
arrange(model.data, AICc)
#save
write.csv(arrange(model.data, AICc), "DFA model selection NOT DETRENDED community all years.csv")

# fit the best model for the 1965-2012 time series
model.list = list(A="zero", m=3, R="diagonal and unequal")

mod = MARSS(all.comm.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# and rotate the loadings
Z.est = coef(mod, type="matrix")$Z
H.inv = varimax(coef(mod, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H.inv)
Z.rot$names <- rownames(all.comm.dat)

# and rotate the trends!
trends.rot = solve(H.inv) %*% mod$states

# DW test for model residuals

# get residuals for each TS
res <- residuals(mod)$residuals
# replace 0 with NA
drop <- res==0
res[drop] <- NA

dw <- dw.p <- NA
for(i in 1:nrow(res)){

  dw.p[i] <- dwtest(res[i,] ~ 1)$p.value
  dw[i] <- dwtest(res[i,] ~ 1)$statistic
}

model.res <- data.frame(names=rownames(res), dw=dw, p=dw.p)
com.mod1.res <- model.res

# try with second-best model
# fit the best model for the 1965-2012 time series
model.list = list(A="zero", m=2, R="diagonal and unequal")

mod2 = MARSS(all.comm.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get residuals for each TS
res <- residuals(mod2)$residuals
# replace 0 with NA
drop <- res==0
res[drop] <- NA

dw <- dw.p <- NA
for(i in 1:nrow(res)){
  
  dw.p[i] <- dwtest(res[i,] ~ 1)$p.value
  dw[i] <- dwtest(res[i,] ~ 1)$statistic
}

model.res <- data.frame(names=rownames(res), dw=dw, p=dw.p)
com.mod2.res <- model.res

####
# third-best
model.list = list(A="zero", m=3, R="equalvarcov")

mod3 = MARSS(all.comm.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get residuals for each TS
res <- residuals(mod3)$residuals
# replace 0 with NA
drop <- res==0
res[drop] <- NA

dw <- dw.p <- NA
for(i in 1:nrow(res)){
  
  dw.p[i] <- dwtest(res[i,] ~ 1)$p.value
  dw[i] <- dwtest(res[i,] ~ 1)$statistic
}

model.res <- data.frame(names=rownames(res), dw=dw, p=dw.p)
com.mod3.res <- model.res
                                         
# now calculate difference between predicted (trend * loading) and observed for each trend!
t1.r <- t2.r <- t3.r <- t4.r <- matrix(nrow=48, ncol=14)
# t.r = matrix(nrow=48, ncol=14)
# pred = matrix(nrow=48, ncol=14)
for(i in 1:14){
  # i <- 1
  t1.r[,i] <- (trends.rot[1,] * Z.rot[i,1]) - all.comm.dat[i,]
  t2.r[,i] <- (trends.rot[2,] * Z.rot[i,2]) - all.comm.dat[i,]
  t3.r[,i] <- (trends.rot[3,] * Z.rot[i,3]) - all.comm.dat[i,]
  #  t4.r[,i] <- (trends.rot[4,] * Z.rot[i,4]) - all.dat[i,]
  # pred[,i] <- (as.matrix(Z.rot[,1:4]) %*% trends.rot)[i,]
  # t.r[,i] <- pred[,i] - all.dat[i,]
}

colnames(t1.r) <- colnames(t2.r) <- colnames(t3.r) <- colnames(t4.r) <- rownames(all.comm.dat)

# colnames(t.r) <- rownames(all.dat)


# calculate AR(1) values for the resulting time series
# using Durbin-Watson test
library(lmtest)
dw.1 <- dw.2 <- dw.3 <- dw.4 <- NA

for(i in 1:14){
  dw.1[i] <- dwtest(t1.r[,i] ~ 1)$p.value
  dw.2[i] <- dwtest(t2.r[,i] ~ 1)$p.value
  dw.3[i] <- dwtest(t3.r[,i] ~ 1)$p.value
  # dw.4[i] <- dwtest(t4.r[,i] ~ 1)$p.value
}

# and get the ar values
ar.1 <- ar.2 <- ar.3 <- NA

for(i in 1:14){
  ar.1[i] <- ar(na.omit(t1.r[,i]), order.max=1, aic=F)$ar
  ar.2[i] <- ar(na.omit(t2.r[,i]), order.max=1, aic=F)$ar
  ar.3[i] <- ar(na.omit(t3.r[,i]), order.max=1, aic=F)$ar
}

comm.all.ar <- data.frame(ar=c(ar.1, ar.2, ar.3))

# combing ar values and plot
comm.65.88.ar$model <- "1965-1988"
comm.89.12.ar$model <- "1989-2012"
comm.all.ar$model <- "1965-2012"

plot.ar <- rbind(comm.65.88.ar, comm.89.12.ar, comm.all.ar)
ggplot(plot.ar, aes(x=ar)) + geom_density(aes(fill = model), alpha=0.4)

# restrict to significant results and combine!
keep <- dw.1 <= 0.05
t1.r <- t1.r[,keep]
t1.r <- gather(as.data.frame(t1.r))
t1.r$year <- 1965:2012
t1.r$trend <- "Trend 1"

keep <- dw.2 <= 0.05
t2.r <- t2.r[,keep]
t2.r <- gather(as.data.frame(t2.r))
t2.r$year <- 1965:2012
t2.r$trend <- "Trend 2"

keep <- dw.3 <= 0.05
t3.r <- t3.r[,keep]
t3.r <- gather(as.data.frame(t3.r))
t3.r$year <- 1965:2012
t3.r$trend <- "Trend 3"

# how many show significant AR(1)?
sum(dw.1<=0.05)
sum(dw.2<=0.05)
sum(dw.3<=0.05)
res.plot <- rbind(t1.r, t2.r, t3.r)

# make a loadings plot!
#Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2,3)])
Z.rot$names <- rownames(all.comm.dat)
Z.rot$key <- rep(c("Trend 1", "Trend 2", "Trend 3"), each=14)
#Z.rot$names <- reorder(Z.rot$names, 1:14)


load.plot <- ggplot(Z.rot, aes(names, value, fill=key)) + geom_bar(stat="identity", position="dodge") +
  theme_gray() + ylab("Loading") + xlab("") + facet_wrap(~key, scales="free_y", nrow=1) +
  scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3], "Trend 3" = cb[4])) +
  geom_hline(yintercept = 0) + theme(legend.position="none") + 
  theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=12)) 

png("community dfa: trends 1-3 predicted-observed NOT DETRENDED.png", 9,6, units="in", res=300)
ggplot(res.plot, aes(year, value, color=trend)) + geom_line() + geom_point() + geom_hline(yintercept = 0, lwd=0.4) +
  facet_wrap(~key, scales="free_y") + theme_gray() + ylab("Trend-specific residual") +
  xlab("") + theme(legend.position = 'top') + theme(legend.title=element_blank())
dev.off()

# and a version with the legend on the inside!

pp <- ggplot(res.plot, aes(year, value, color=trend)) + geom_line() + geom_point() +
  facet_wrap(~key, scales="free_y") + theme_gray() + ylab("Trend-specific residual") + 
  xlab("") + theme(legend.title=element_blank()) + scale_color_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3], "Trend 3" = cb[4]))

png("non-detrended community dfa: trends 1-3 residuals.png", 8,6, units="in", res=300)
reposition_legend(pp, 'left', panel = 'panel-4-2')
dev.off()

other.plot <- reposition_legend(pp, 'left', panel = 'panel-4-2')

png("non-detrended community dfa: loading and residuals.png", 8,10, units="in", res=300)
ggarrange(load.plot,  other.plot, labels = c("a)", "b)"), nrow=2, heights = c(0.9,1))
dev.off()



######################################################################
## and, now figure out the trend difference in slope and intercept!
######################################################################
# re-rotate the loadings, which were messed up in the plot code above!
# and rotate the loadings
Z.est = coef(mod, type="matrix")$Z
H.inv = varimax(coef(mod, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H.inv)
Z.rot$names <- rownames(all.comm.dat)

# re-calculate residuals...
t1.r <- t2.r <- t3.r <- matrix(nrow=48, ncol=14)
# t.r = matrix(nrow=48, ncol=14)
# pred = matrix(nrow=48, ncol=14)
for(i in 1:14){
#   i <- 1
  t1.r[,i] <- (trends.rot[1,] * Z.rot[i,1]) - all.comm.dat[i,]
  t2.r[,i] <- (trends.rot[2,] * Z.rot[i,2]) - all.comm.dat[i,]
  t3.r[,i] <- (trends.rot[3,] * Z.rot[i,3]) - all.comm.dat[i,]
  #  t4.r[,i] <- (trends.rot[4,] * Z.rot[i,4]) - all.dat[i,]
  # pred[,i] <- (as.matrix(Z.rot[,1:4]) %*% trends.rot)[i,]
  # t.r[,i] <- pred[,i] - all.dat[i,]
}

colnames(t1.r) <- colnames(t2.r) <- colnames(t3.r) <- colnames(t4.r) <- rownames(all.comm.dat)

# restrict to significant results
keep <- dw.1 <= 0.05
t1.r <- t1.r[,keep]

keep <- dw.2 <= 0.05
t2.r <- t2.r[,keep]

keep <- dw.3 <= 0.05
t3.r <- t3.r[,keep]

rownames(t1.r) <- rownames(t2.r) <- rownames(t3.r) <- 1965:2012
start <- 1980:1998

tr1.aicc <- tr2.aicc <- matrix(nrow=19, ncol=ncol(t1.r))

for(i in 1:19){
 # i <- 1
  
  for(j in 1:ncol(t1.r)){
 # j <- 1  
  temp <- data.frame(res=t1.r[,j], year=1965:2012, era=1)  
  temp$era[temp$year >= start[i]] <- 2  
  mod <- lm(res ~ year*era, data=temp)
  tr1.aicc[i,j] <- AICc(mod)
  
  temp <- data.frame(res=t2.r[,j], year=1965:2012, era=1)  
  temp$era[temp$year >= start[i]] <- 2  
  mod <- lm(res ~ year*era, data=temp)
  tr2.aicc[i,j] <- AICc(mod) 
  
  
  }
}

tr3.aicc <- matrix(nrow=19, ncol=ncol(t3.r))

for(i in 1:19){
  # i <- 1
  
  for(j in 1:ncol(t3.r)){
temp <- data.frame(res=t3.r[,j], year=1965:2012, era=1)  
temp$era[temp$year >= start[i]] <- 2  
mod <- lm(res ~ year*era, data=temp)
tr3.aicc[i,j] <- AICc(mod)
}}

plot(start-0.5, rowMeans(tr1.aicc), type="l")
plot(start-0.5, rowMeans(tr2.aicc), type="l")
plot(start-0.5, rowMeans(tr3.aicc), type="l")

tr1.daicc <- rowMeans(tr1.aicc)-min(rowMeans(tr1.aicc))
tr2.daicc <- rowMeans(tr2.aicc)-min(rowMeans(tr2.aicc))
tr3.daicc <- rowMeans(tr3.aicc)-min(rowMeans(tr3.aicc))

plot.all <- data.frame(aic=c(tr1.daicc, tr2.daicc, tr3.daicc), year=start-0.5, 
                       trend=rep(c("Trend 1", "Trend 2", "Trend 3"), each=19))

# get AICc for null model (stationary residuals) to compare
null.1 <- NA

for(j in 1:ncol(t1.r)){
  temp <- data.frame(res=t1.r[,j])  
  mod <- lm(res ~ 1, data=temp)
  null.1[j] <- AICc(mod)
}
names(null.1) <- colnames(t1.r)
##
null.2 <- NA

for(j in 1:ncol(t2.r)){
  temp <- data.frame(res=t2.r[,j])  
  mod <- lm(res ~ 1, data=temp)
  null.2[j] <- AICc(mod)
}
names(null.2) <- colnames(t2.r)
##
null.3 <- NA

for(j in 1:ncol(t3.r)){
  temp <- data.frame(res=t3.r[,j])  
  mod <- lm(res ~ 1, data=temp)
  null.3[j] <- AICc(mod)
}
names(null.3) <- colnames(t3.r)

# compare with minimum AICc from nontationary model
nst1.best <- nst2.best <- nst3.best <- NA

for(i in 1:length(null.1)){

  ifelse(null.1[i] > mean(tr1.aicc[,i]), nst1.best[i] <- 1,  nst1.best[i] <- 1)
         }
sum(nst1.best); length(nst1.best)
##
for(i in 1:length(null.2)){
  
  ifelse(null.2[i] > mean(tr2.aicc[,i]), nst2.best[i] <- 1,  nst2.best[i] <- 1)
}
sum(nst2.best); length(nst2.best)
##
for(i in 1:length(null.3)){
ifelse(null.3[i] > mean(tr3.aicc[,i]), nst3.best[i] <- 1,  nst3.best[i] <- 1)
}
sum(nst3.best); length(nst3.best)

pdf("aicc for model residuals.pdf", 4,3.5)
ggplot(plot.all, aes(year, aic, color=trend)) + geom_line() + theme_gray() + 
  theme(legend.title = element_blank(), legend.position = c(0.2,0.8)) + xlab("") + ylab("delta-AICc")
dev.off()

# and combine with the loadings and residuals!
aic.plot <- ggplot(plot.all, aes(year, aic, color=trend)) + geom_line() + theme_gray() + 
  scale_color_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3], "Trend 3" = cb[4])) +
  theme(legend.title = element_blank(), legend.position = c(0.3,0.75)) + xlab("") + ylab("delta-AICc")

half.plot <- ggarrange( aic.plot, labels = c("c)"), ncol=3)

png("non-detrended community dfa: loading residuals and aic.png", 8,12, units="in", res=300)
ggarrange(load.plot,  other.plot, half.plot, labels = c("a)", "b)"), nrow=3, heights = c(0.9,1,0.5), widths = c(1,1,0.7))
dev.off()

# and a restricted version with no loadings
half.plot <- ggarrange( aic.plot, labels = c("b)"), ncol=3)
png("non-detrended community dfa: residuals and aic NO loadings.png", 8,8, units="in", res=300)
ggarrange(other.plot, half.plot, labels = "a)", nrow=2, heights = c(1.1,0.5), widths = c(1,0.7))
dev.off()
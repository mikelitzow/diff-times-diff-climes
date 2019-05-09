# this is a streamlined version of DFA analysis that was used in the submitted version of the manuscript

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
library(lmtest)
library(cowplot)

# DFA code for paper                                                                     
#
# this code does two things:
#
# first, fit DFA models to biology and environmental data from
# the two eras - before and after 1988/89
# 
# second, fit DFA models to the entire data set and calculate rolling window correlations
# between shared trends and individual time series
#########################################################################################
# # begin with environmental data
# # load data...
# dat <- read.csv("/Users/MikeLitzow 1/Documents/R/time and climes/clim vars for ordination.csv", row.names=1)
# 
# # put years into row names
# rownames(dat) <- dat[,1]
# dat <- dat[,-1]
# 
# # dropping AL and NPI as I think we will leave out these dynamics to keep the paper simple
# # the relevant atmosphere-ocean dynamics have been addequately addressed in PRSB ms. for our purposes!
# dat <- dat[,-c(1,2)]
# 
# head(dat)
# 
# # remove some extra  time series
# dat <- dat[,-c(5,6,8,9,10,11,14,17)]
# 
# 
# # change names to plot-friendly labels
# colnames(dat) <- c("PDO", "NPGO", "SLP gradient", "Freshwater", "Wind stress", "Downwelling", "SST", "Advection", "SSH")
# 
# write.csv(dat, "GOA environmental data.csv")

# load environmental data
dat <- read.csv("GOA environmental data.csv", row.names = 1)

# change names to plot-friendly labels
colnames(dat)[c(3,5)] <- c("SLP gradient", "Wind stress")

# plot TS for SI
do.dat <- as.data.frame(scale(dat)) # scale to plot on 1 axis
plot.dat <- gather(do.dat)
plot.dat$year <- 1950:2012

# make an order key for plotting
plot.dat$order <- rep(c(9,8,4,3,7,2,5,1,4), each=nrow(dat))

# and reorder labels
plot.dat$key <- reorder(plot.dat$key, plot.dat$order)

# save as SI Fig
# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf("SI - climate time series.pdf", 8,6)
ggplot(plot.dat, aes(x=year, y=value)) + geom_bar(position="dodge", stat="identity", fill=cb[2]) + facet_wrap(~key) +
  ylab("Standard anomaly") + xlab("") + theme_bw() + geom_hline(yintercept = 0)
dev.off()

# remove PDO/NPGO as they are not included in the analysis of regional environmental variability
dat <- dat[,-c(1,2)]

# check
head(dat)

# save the full data set for later use...
all.clim.dat <- t(as.matrix(dat))

# now fit DFA models
# first the early era
e.cli.dat = as.matrix(dat[rownames(dat) %in% 1950:1988,])

# and transpose
e.cli.dat <- t(e.cli.dat)

# now fit DFA models with 1-3 trends and different error structures and compare

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# fit models & store results
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

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data

# 2 trend equal var covar is the best model....
# save the model output table for SI
write.csv(model.data, "climate dfa table 1950-1988.csv")
early.model.data <- model.data

############
# and the late era!

l.cli.dat = as.matrix(dat[rownames(dat) %in% 1989:2012,])

# and transpose
l.cli.dat <- t(l.cli.dat)
model.data = data.frame()

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
write.csv(model.data, "climate dfa table 1989-2012.csv")

############
# now fit best model for each era

model.list = list(A="zero", m=2, R="equalvarcov") # best model for early era
mod.early = MARSS(e.cli.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# and rotate the loadings
Z.est = coef(mod.early, type="matrix")$Z
H.inv = varimax(coef(mod.early, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H.inv)

# reverse trend 2 to plot
Z.rot[,2] <- -Z.rot[,2]

Z.rot$names <- rownames(e.cli.dat)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(e.cli.dat)
Z.rot$plot.names <- reorder(Z.rot$names, 1:14)

# and the late era
# equal var covar one trend is the best for this era!
model.list = list(A="zero", m=1, R="equalvarcov")
mod.late = MARSS(l.cli.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...
modCI.late <- MARSSparamCIs(mod.late)

plot.CI.late <- data.frame(names=rownames(l.cli.dat), mean=modCI.late$par$Z, upCI=modCI.late$par.upCI$Z,
                           lowCI=modCI.late$par.lowCI$Z)

plot.CI.late <- arrange(plot.CI.late, mean)
plot.CI.late$names.order <- reorder(plot.CI.late$names, plot.CI.late$mean)
dodge <- position_dodge(width=0.9)

# combine into one plot
Z.rot$key <- rep(c("Trend 1", "Trend 2"), each=7)
plot.CI.late$names.order <- reorder(plot.CI.late$names, plot.CI.late$mean)


early.env <- ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar(stat="identity", position="dodge") +
  theme_bw() + ylab("Loading") + xlab("") + ggtitle("1950-1988 environment") + 
  scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3])) +
  theme(legend.position = c(0.8,0.2), legend.title=element_blank()) + geom_hline(yintercept = 0) +
  theme(axis.text.x  = element_text(angle=45, hjust=1, size=12)) + ylim(-0.6, 0.8)

# reorder Z plot
plot.CI.late$plot.names <- reorder(plot.CI.late$names, c(1,3,5,2,6,7,4))

late.env <- ggplot(plot.CI.late, aes(x=plot.names, y=mean)) + geom_bar(position="dodge", stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme_bw() + 
  theme(axis.text.x  = element_text(angle=45, hjust=1,  size=12), axis.title.y = element_blank()) +
  geom_hline(yintercept = 0) + ggtitle("1989-2012 environment") + ylim(-0.6, 0.8)

# # combine
# png("loadings plot non-detrended environmental dfa both eras.png", 6, 4, units="in", res=300)
# ggarrange(early.plot,  late.plot, labels = c("a)", "b)"),  widths=c(1.3,1), ncol=2)
# dev.off()


##############
# now find the best model for entire environmental time series


# cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# object for model data
model.data = data.frame()

# object for residual autocorrelation results
res.ar <- matrix(nrow=1, ncol=10)

# fit models
for(R in levels.R) {
  for(m in 1:3) {
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(all.clim.dat, model=dfa.model, control=cntl.list,
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
    
    # drop the 0s - these should be NAs!
    
    drop <- res==0
    res[drop] <- NA
    
    for(ii in 1:nrow(res)){
      
      dw.p[ii] <- dwtest(res[ii,] ~ 1)$p.value
      
    }
    
    # pad to the correct length for 1- and 2-trend models
    if(length(dw.p)==8) {dw.p <- c(dw.p, NA, NA)}
    if(length(dw.p)==9) {dw.p <- c(dw.p, NA)}
    
    res.ar <- rbind(res.ar, dw.p)
    
  } # end m loop
} # end R loop

res.ar
colnames(res.ar) <- rownames(res) # make sure this case of "res" is a full case, i.e., three shared trends
res.ar <- res.ar[2:nrow(res.ar),] # drop first row of NAs

model.res.table <- as.data.frame(res.ar)
model.res.table$R <- model.data$R
model.res.table$m <- model.data$m

write.csv(model.res.table, "residual ar1 all environmental models 1965-2012.csv")

###
arrange(model.data, AICc)

# save model table
write.csv(arrange(model.data, AICc), "DFA model selection environment all years.csv")

#########
# fit the best model to the entire time series
# and calculate the moving correlations for each TS!

cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
model.list = list(A="zero", m=1, R="unconstrained")
mod <- MARSS(all.clim.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

trend1 <- as.data.frame(matrix(nrow=24, ncol=7))
colnames(trend1) <- rownames(all.clim.dat)

for(i in 1950:1988){
  # i <- 1950
  temp <- all.clim.dat[,colnames(all.clim.dat) %in% i:(i+24)]
  
  for(ii in 1:7){ # loop through each time series/variable
    
    trend1[(i-1949),ii] <- cor(temp[ii,], mod$states[1,(i-1949):(i-1949+24)], use="p")
    
  } }

# now  plot

plot1 <- gather(trend1)
plot1$year <- 1962:2000

# plot2 <- gather(trend2)
# plot2$year <- 1977:2000
# plot2$trend <- "Trend 2"
# 
# plot3 <- gather(trend3)
# plot3$year <- 1977:2000
# plot3$trend <- "Trend 3"

# # set colors
# cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# png("environmental correlations with trend for best DFA model entire time series.png", 8,3, units="in", res=300)
# ggplot(plot1, aes(year, value)) +
#   facet_wrap(~key, nrow=2, ncol=4) + geom_hline(yintercept = 0, color="dark grey") + theme_grey() +
#   geom_line(color=cb[2]) + 
#   xlab("") + ylab("Correlation") + theme(axis.title.x = element_blank())
# dev.off()

# save for combined plot!
env.corr.plot <- ggplot(plot1, aes(year, value)) +
  facet_wrap(~key, nrow=2, ncol=4) + geom_hline(yintercept = 0) + theme_bw() +
  geom_line(color=cb[2]) + 
  xlab("") + ylab("Correlation") + theme(axis.title.x = element_blank()) + 
  xlim(1962,2000) + 
  ggtitle("Environment")

####
# and fit the second-best model, which has independent residuals!
model.list = list(A="zero", m=1, R="equalvarcov")
mod <- MARSS(all.clim.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

trend1 <- as.data.frame(matrix(nrow=24, ncol=7))
colnames(trend1) <- colnames(trend2) <- colnames(trend3) <- colnames(dat)

for(i in 1950:1988){
  # i <- 1950
  temp <- t(dat[rownames(dat) %in% i:(i+24),])
  
  for(ii in 1:7){ # loop through each time series/variable
    
    trend1[(i-1949),ii] <- cor(temp[ii,], mod$states[1,(i-1949):(i-1949+24)], use="p")
    
  } }

# now  plot

plot1 <- gather(trend1)
plot1$year <- 1962:2000

# png("environmental correlations with trend for SECOND-best DFA model entire time series.png", 8,3, units="in", res=300)
# ggplot(plot1, aes(year, value)) +
#   facet_wrap(~key, nrow=2, ncol=4, scales="free_y") + theme_grey() +
#   geom_line(color=cb[2]) + 
#   xlab("") + ylab("Correlation") + theme(axis.title.x = element_blank())
# dev.off()

# save for SI

pdf("SI - environment corrs 2nd best dfa model.pdf", 8,3)
ggplot(plot1, aes(year, value)) +
  facet_wrap(~key, nrow=2, ncol=4, scales="free_y") + theme_bw() +
  geom_line(color=cb[2]) + 
  xlab("") + ylab("Correlation") + theme(axis.title.x = element_blank())
dev.off()

################################
# now...the community analysis #
################################
# 
# dat <- read.csv("salmon and non-salmon biology mar 28.csv", row.names = 1)
# 
# # drop GOAPOP
# keep <- colnames(dat)!="GOAPOP"
# dat <- dat[,keep] 
# 
# # check
# head(dat)
# 
# # limit to 1965:2012
# dat <- dat[rownames(dat) <= 2012,]
# 
# write.csv(dat, "GOA community data.csv")

dat <- read.csv("GOA community data.csv", row.names=1)

# change names to plot-friendly labels
colnames(dat) <- c("Sablefish", "Walleye pollock", "Arrowtooth flounder", "Pacific herring", "Tanner crab", "Sidestripe shrimp", "Capelin", "Northern shrimp", "Pacific cod", "Chum salmon", "Pink salmon",
                   "Coho salmon", "Sockeye salmon", "Northern rockfish")

# plot time series for SI
do.dat <- as.data.frame(scale(dat)) # scale to plot on 1 axis
plot.dat <- gather(do.dat)
plot.dat$year <- 1965:2012

# save as SI Fig

pdf("SI - community time series.pdf", 8,6)
ggplot(plot.dat, aes(x=year, y=value)) + geom_bar(position="dodge", stat="identity", fill=cb[2]) +
  facet_wrap(~key) +
  ylab("Standard anomaly") + xlab("") + geom_hline(yintercept = 0) + theme_bw() 
dev.off()

# fit to the early era
e.comm.dat = as.matrix(dat[rownames(dat) %in% 1965:1988,])

# and transpose
e.comm.dat <- t(e.comm.dat)

# fit models
model.data = data.frame()

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

# save results
write.csv(arrange(model.data, AICc), "DFA model selection community 65-88.csv")

early.model.data <- model.data

# and fit to the late era
l.comm.dat = as.matrix(dat[rownames(dat) %in% 1989:2012,])

# and transpose
l.comm.dat <- t(l.comm.dat)

model.data = data.frame()

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

# save results
write.csv(arrange(model.data, AICc), "DFA model selection community 89-12.csv")

###########
# fit best model in each era
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
# ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar(stat="identity", position="dodge") +
#   theme_gray() 

# # and rotate the trends!
# trends.rot = solve(H.inv) %*% mod.early$states
# 
# # plot these trends!
# plot.trends <- as.data.frame(t(trends.rot)) %>%
#   gather()
# plot.trends$year <- 1965:1988
# 
# ggplot(plot.trends, aes(year, value, color=key)) + geom_line() + theme_gray()

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

# ggplot(plot.CI.late, aes(x=names.order, y=mean)) + geom_bar(position="dodge", stat="identity") +
#   geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
#   ylab("Loading") + xlab("") + theme(axis.text.x = element_text(angle = 90, vjust=0, hjust=1)) +
#   geom_hline(yintercept = 0) + theme_gray()

# combine into one plot
Z.rot$key <- rep(c("Trend 1", "Trend 2"), each=14)
#plot.CI.late$names.order <- reorder(plot.CI.late$names, plot.CI.late$mean)

early.comm <- ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar( position="dodge",stat="identity") +
  theme_bw() + ylab("Loading") + xlab("") + ggtitle("1965-1988 community") + 
  scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3])) +
  theme(legend.position = c(0.12,0.85), legend.title=element_blank()) + geom_hline(yintercept = 0) +
  theme(axis.text.x  = element_text(angle=45, hjust=1, size=12)) +
  ylim(-0.47,0.56)

late.comm <- ggplot(plot.CI.late, aes(x=names.order, y=mean)) + 
  geom_bar(position="dodge", stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("") + xlab("") + theme_bw() + theme(axis.text.x  = element_text(angle=45, hjust=1, size=12)) +
  geom_hline(yintercept = 0) + ggtitle("1989-2012 community") + 
  ylim(-0.47,0.56)

# # combine
# png("loadings plot non-detrended community dfa both eras.png", 8, 5, units="in", res=300)
# ggarrange(early.plot,  late.plot, labels = c("c)", "d)"),  widths=c(1.2,1), ncol=2)
# dev.off()

# combine era-specific environment and biology loadings
# into a single plot

# start by making a blank object for spacing the top row
blank <- ggplot() + theme_void()

top <- plot_grid(blank, early.env, late.env, blank,
                 labels = c(NA, 'a)', 'b)', NA), rel_widths = c(0.15,1.1, 0.7, 0.25), nrow=1)

bottom <- plot_grid(early.comm, late.comm,
                 labels = c('c)', 'd)'), rel_widths = c(1.2, 1), nrow=1)

pdf("environment and biology loadings by era.pdf", 9,9)
plot_grid(top, bottom, nrow=2, rel_heights = c(0.9, 1))
dev.off()


# get correlations between correlations in the two eras
# make a data frame to compbine all loadings
cor.df <- data.frame(names=Z.rot$names[1:14], load1=Z.rot$value[Z.rot$key=="Trend 1"], load2=Z.rot$value[Z.rot$key=="Trend 2"])

# join with late era model loadings
cor.df=left_join(cor.df, plot.CI.late)

cor.test(cor.df$load1, cor.df$mean)
# Pearson's product-moment correlation
# 
# data:  cor.df$load1 and cor.df$mean
# t = -0.71585, df = 12, p-value = 0.4878
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.6618824  0.3676876
# sample estimates:
# cor 
# -0.2023722 

cor.test(cor.df$load2, cor.df$mean)
# Pearson's product-moment correlation
# 
# data:  cor.df$load2 and cor.df$mean
# t = -1.2337, df = 12, p-value = 0.2409
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.7352094  0.2373185
# sample estimates:
# cor 
# -0.3355068

###############################################
# now find the best model for entire community time series
all.comm.dat = dat

# and transpose
all.comm.dat <- t(all.comm.dat)

cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# fit models and calculate residual autocorrelation
model.data = data.frame()
res.ar <- matrix(nrow=1, ncol=17)

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


###
arrange(model.data, AICc)
#save
write.csv(arrange(model.data, AICc), "DFA model selection community all years.csv")

#########
# fit the best model to the entire time series
# and calculate the moving correlations for each TS!

model.list = list(A="zero", m=3, R="diagonal and unequal")
mod <- MARSS(all.comm.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# and rotate the loadings
Z.est = coef(mod, type="matrix")$Z
H.inv = varimax(coef(mod, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H.inv)
Z.rot$names <- colnames(dat)

# and rotate the trends!
trends.rot = solve(H.inv) %*% mod$states

trend1 <- trend2 <- trend3 <- as.data.frame(matrix(nrow=24, ncol=14))
colnames(trend1) <- colnames(trend2) <- colnames(trend3) <- colnames(dat)

for(i in 1965:1988){
 # i <- 1965
  temp <- t(dat[rownames(dat) %in% i:(i+24),])
  
 for(ii in 1:14){ # loop through each population
   
   trend1[(i-1964),ii] <- cor(temp[ii,], trends.rot[1,(i-1964):(i-1964+24)], use="p")
   trend2[(i-1964),ii] <- cor(temp[ii,], trends.rot[2,(i-1964):(i-1964+24)], use="p")
   trend3[(i-1964),ii] <- cor(temp[ii,], trends.rot[3,(i-1964):(i-1964+24)], use="p")
 } }
  
# now combine and plot

plot1 <- gather(trend1)
plot1$year <- 1977:2000
plot1$trend <- "Trend 1"

plot2 <- gather(trend2)
plot2$year <- 1977:2000
plot2$trend <- "Trend 2"

plot3 <- gather(trend3)
plot3$year <- 1977:2000
plot3$trend <- "Trend 3"
plot <- rbind(plot1, plot2, plot3)


pp <- ggplot(plot, aes(year, value, color=trend)) + geom_line() + 
  facet_wrap(~key) + geom_hline(yintercept = 0) + theme_bw() +
  xlab("") + ylab("Correlation") + theme(legend.title=element_blank(), axis.title.x = element_blank())+ 
  scale_color_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3], "Trend 3" = cb[4])) +
  xlim(1962,2000) +
  ggtitle("Community")#+

comm.corr.plot <- reposition_legend(pp, 'center', panel = 'panel-4-3')

pdf("combined environment and community correlations with trends 1-3.pdf", 8,9)
ggarrange(env.corr.plot,  comm.corr.plot, labels = c("a)", "b)"), nrow=2, heights=c(1.1,2))
dev.off()


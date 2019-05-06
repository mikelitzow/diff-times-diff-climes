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

# begin with environmental data
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
  theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=12)) + ylim(-0.6, 0.8)

# reorder Z plot
plot.CI.late$plot.names <- reorder(plot.CI.late$names, c(1,3,5,2,6,7,4))

late.plot <- ggplot(plot.CI.late, aes(x=plot.names, y=mean)) + geom_bar(position="dodge", stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme_gray() + theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=12)) +
  geom_hline(yintercept = 0) + ggtitle("1989-2012") + ylim(-0.6, 0.8) + ylab("")

# combine
pdf("loadings plot detrended environmental dfa both eras.pdf", 6, 4)
ggarrange(early.plot,  late.plot, labels = c("a)", "b)"),  widths=c(1.3,1), ncol=2)
dev.off()


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
write.csv(arrange(model.data, AICc), "DFA model selection community 65-88.csv")

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
write.csv(arrange(model.data, AICc), "DFA model selection community 89-12.csv")

###############################################
# fit best model for each era
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

model.list = list(A="zero", m=2, R="diagonal and unequal") # best model for early era
mod.early = MARSS(e.comm.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# need to rotate!
# and rotate the loadings
Z.est = coef(mod.early, type="matrix")$Z
H.inv = varimax(coef(mod.early, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H.inv)

Z.rot$names <- rownames(e.comm.dat)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(e.comm.dat)
Z.rot$plot.names <- reorder(Z.rot$names, 28:1)
ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar(stat="identity", position="dodge") +
  theme_gray() 

plot.early <- Z.rot
#################################
# and the late model

# equal var covar one trend is the best for this era!
model.list = list(A="zero", m=1, R="diagonal and equal")
mod.late = MARSS(l.comm.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...
modCI.late <- MARSSparamCIs(mod.late)
modCI.late

# combine into one plot
Z.rot$key <- rep(c("Trend 1", "Trend 2"), each=7)
plot.CI.late$names.order <- reorder(plot.CI.late$names, plot.CI.late$mean)

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

early.plot <- ggplot(Z.rot, aes(plot.names, value, fill=key)) + geom_bar(stat="identity", position="dodge") +
  theme_gray() + ylab("Loading") + xlab("") + ggtitle("1950-1988") + 
  scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3])) +
  theme(legend.position = c(0.8,0.2), legend.title=element_blank()) + geom_hline(yintercept = 0) +
  theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=12)) + ylim(-0.6, 0.8)

# reorder Z plot
plot.CI.late$plot.names <- reorder(plot.CI.late$names, c(1,3,5,2,6,7,4))

late.plot <- ggplot(plot.CI.late, aes(x=plot.names, y=mean)) + geom_bar(position="dodge", stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme_gray() + theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=12)) +
  geom_hline(yintercept = 0) + ggtitle("1989-2012") + ylim(-0.6, 0.8) + ylab("")

# combine
pdf("loadings plot detrended environmental dfa both eras.pdf", 6, 4)
ggarrange(early.plot,  late.plot, labels = c("a)", "b)"),  widths=c(1.3,1), ncol=2)
dev.off()


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


###
arrange(model.data, AICc)
#save
write.csv(arrange(model.data, AICc), "DFA model selection community all years.csv")

#########
# fit the best model to the entire time series
# and calculate the moving correlations for each TS!

cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
model.list = list(A="zero", m=3, R="diagonal and unequal")
mod <- MARSS(t(dat), model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

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

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pp <- ggplot(plot, aes(year, value, color=trend)) + geom_line() + 
  facet_wrap(~key) + geom_hline(yintercept = 0, color="dark grey") + theme_gray() +
  xlab("") + ylab("Correlation") + theme(legend.title=element_blank())+ 
  scale_color_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3], "Trend 3" = cb[4]))

png("community correlations with trends 1-3.png", 8,6, units="in", res=300)
reposition_legend(pp, 'center', panel = 'panel-4-3')
dev.off()

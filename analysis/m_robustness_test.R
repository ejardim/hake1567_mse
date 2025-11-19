#==================================================================
# script for M analysis
# EJ
# 20251027
#==================================================================


#------------------------------------------------------------------
# load & set up
#------------------------------------------------------------------

library(FLa4a)
library(ggplotFL)
library(reshape2)
library(knitr)
# stock
load("../data/HKE_1_5_6_7_stk_input_assess.Rdata")
ages  <- as.numeric(dimnames(m(hke.stk))$age)
years <- as.numeric(dimnames(m(hke.stk))$year)
rng <- range(hke.stk)
# weights at age
wa <- yearMeans(stock.wt(hke.stk))
# index
load("../data/FLIndices_0.RData") # medits data fixed
hke.idx <- trim(idx[[1]], age=0:4, year=years)
# number of cores
cores <- 5

#------------------------------------------------------------------
# growth model and parameters
#------------------------------------------------------------------

t0 <- -0.005
k <- 0.178
linf <- 110
a <- 0.00677
b <- 3.0325
lmin <- 5
l50 <- 29
maxage <- t0-1/k*log(1-(linf-0.5)/linf)
wa <- yearMeans(stock.wt(hke.stk))
#VBGF function Bertalanffy Growth Function (VBGF):
hke1567.vb <- a4aGr(
  grMod=~linf*(1-exp(-k*(t-t0))),
  grInvMod=~t0-1/k*log(1-len/linf),
  params=FLPar(linf=linf, k=k, t0=t0, units=c("cm","year-1","year"))
)

#------------------------------------------------------------------
# M models
#------------------------------------------------------------------

# current model
m.c <- m(hke.stk)

# age-dependent version of the Chen & Watanabe model, Chen, S., & Watanabe, S. (1989)
mChen_sim <- FLModelSim(model=~k / (1 - exp(-k * (age - t0))),
  params=FLPar(k=k, t0 = t0))
m.cw <- m.c
m.cw[] <- predict(mChen_sim, age=ages+0.5)

# constant M, Then, A. Y., Hoenig, J. M., Hall, N. G., & Hewitt, D. A. (2015)
mThen_sim <- FLModelSim(model=~4.899 * max_age^(-0.916),
  params=FLPar(max_age=maxage))
a4aM_Then <- a4aM(level = mThen_sim)
range(a4aM_Then, c("min","max")) <- rng[c("min", "max")]
range(a4aM_Then, c("minyear","maxyear")) <- rng[c("minyear", "maxyear")]
m.t <- m(a4aM_Then)

# weight based Lorenz
mLoz_sim <- FLModelSim(model=~3*wt^(-0.288))
m.l <- m.c
m.l[] <- predict(mLoz_sim, wt=wa*1000)

# length based Brodziak Brodziak, J., Ianelli, J., Lorenzen, K., & Method, R. D. (2011)
mBrod_sim <- FLModelSim(model=~k*l50/len, params=FLPar(l50=l50, k=k))
mBrod <- FLQuant(dimnames=list(len=lmin:linf, year=rng["minyear"]:rng["maxyear"]))
mBrod[] <- predict(mBrod_sim, len=lmin:linf+0.5)
mBrod <- l2a(mBrod, hke1567.vb, stat="mean")
m.b <- mBrod[ac(ages)]

# plot
df0 <- FLQuants(list("Chen Whatanabe"=m.cw, Brodziak=m.b, Current=m.c, Lorenz=m.l, Then=m.t))
df0 <- lapply(df0, "[", j=1)
df0 <- do.call(cbind, df0)
df0 <- melt(df0, c("age", "model"))
df0$age <- df0$age-0.5
png("../paper/mmodels.png", 800, 800)
xyplot(value ~ age, groups = model, data = df0, type = "b",
  ylab = "M", auto.key = list(space = "top", columns = 5),
  par.settings = list(
    superpose.symbol = list(pch = 19),
    superpose.line = list(lwd = 2)
  )
)
dev.off()

#------------------------------------------------------------------
# Set up a4a stock assessment model #------------------------------------------------------------------

fmod <- ~s(age, k = 4) + s(year, k = 8) + te(age, year, k = c(3, 10))
qmod <- list(~I(1/(1 + exp(-age))))

## mcmc set up
# random seed
rs <- 123
# number of iters
it <- 500
# thinning
mcs <- 1000
# burnin
burnin <- 1000
# control file
mc <- SCAMCMC(mcmc=mcs*(it+burnin), mcsave=mcs, mcseed=rs)

#------------------------------------------------------------------
# OM = current #------------------------------------------------------------------

## Operating model
om.c <- hke.stk
#fit <- sca(om.c, hke.idx, fmodel=fmod, qmodel=qmod, fit="MCMC", mc=mc)
# remove mcscale burnin period
#fit <- burnin(fit, burnin)
fit <- simulate(sca(om.c, hke.idx, fmodel=fmod, qmodel=qmod), it, seed=rs)

# update stock
om.c <- om.c+fit

## OEM
stk.oem <- om.c
idx.oem <- hke.idx
index(idx.oem) <- index(fit)[[1]]
idx.oem <- FLIndices(idx=idx.oem)

## MP

### current
stk.c <- stk.oem
m(stk.c)[] <- m.c

### chen watanabe
stk.cw <- stk.oem
m(stk.cw)[] <- m.cw

### Then
stk.t <- stk.oem
m(stk.t)[] <- m.t

### Lorenz
stk.l <- stk.oem
m(stk.l)[] <- m.l

### Brodziak
stk.b <- stk.oem
m(stk.b)[] <- m.b

### call sca in parallel
stks <- FLStocks(c=stk.c, cw=stk.cw, t=stk.t, l=stk.l, b=stk.b)
idxs <- list(c=idx.oem, cw=idx.oem, t=idx.oem, l=idx.oem, b=idx.oem)
fits <- scas(stks, idxs, fmodel=list(fmod), qmodel=list(qmod), fit="MP", workers=cores)

### output (add om)
output.c <- stks + fits
output.c <- FLStocks(c(output.c, FLStocks(om=om.c)))

#------------------------------------------------------------------
# OM = chen watanabe #------------------------------------------------------------------

## Operating model
om.cw <- hke.stk
m(om.cw)[] <- m.cw
#fit <- sca(om.cw, hke.idx, fmodel=fmod, qmodel=qmod, fit="MCMC", mc=mc)
# remove mcscale burnin period
#fit <- burnin(fit, burnin)
fit <- simulate(sca(om.cw, hke.idx, fmodel=fmod, qmodel=qmod), it, seed=rs)

om.cw <- om.cw+fit

## OEM
stk.oem <- om.cw
idx.oem <- hke.idx
index(idx.oem) <- index(fit)[[1]]
idx.oem <- FLIndices(idx=idx.oem)

## MP

### current
stk.c <- stk.oem
m(stk.c)[] <- m.c

### chen watanabe
stk.cw <- stk.oem
m(stk.cw)[] <- m.cw

### Then
stk.t <- stk.oem
m(stk.t)[] <- m.t

### Lorenz
stk.l <- stk.oem
m(stk.l)[] <- m.l

### Brodziak
stk.b <- stk.oem
m(stk.b)[] <- m.b

### call sca in parallel
stks <- FLStocks(c=stk.c, cw=stk.cw, t=stk.t, l=stk.l, b=stk.b)
idxs <- list(c=idx.oem, cw=idx.oem, t=idx.oem, l=idx.oem, b=idx.oem)
fits <- scas(stks, idxs, fmodel=list(fmod), qmodel=list(qmod), fit="MP", workers=cores)

### output
output.cw <- stks + fits
output.cw <- FLStocks(c(output.cw, FLStocks(om=om.cw)))

#------------------------------------------------------------------
# OM = then #------------------------------------------------------------------

## Operating model
om.t <- hke.stk
m(om.t)[] <- m.t
#fit <- sca(om.t, hke.idx, fmodel=fmod, qmodel=qmod, fit="MCMC", mc=mc)
# remove mcscale burnin period
#fit <- burnin(fit, burnin)
fit <- simulate(sca(om.t, hke.idx, fmodel=fmod, qmodel=qmod), it, seed=rs)

om.t <- om.t+fit

## OEM
stk.oem <- om.t
idx.oem <- hke.idx
index(idx.oem) <- index(fit)[[1]]
idx.oem <- FLIndices(idx=idx.oem)

## MP

### current
stk.c <- stk.oem
m(stk.c)[] <- m.c

### chen watanabe
stk.cw <- stk.oem
m(stk.cw)[] <- m.cw

### Then
stk.t <- stk.oem
m(stk.t)[] <- m.t

### Lorenz
stk.l <- stk.oem
m(stk.l)[] <- m.l

### Brodziak
stk.b <- stk.oem
m(stk.b)[] <- m.b

### call sca in parallel
stks <- FLStocks(c=stk.c, cw=stk.cw, t=stk.t, l=stk.l, b=stk.b)
idxs <- list(c=idx.oem, cw=idx.oem, t=idx.oem, l=idx.oem, b=idx.oem)
fits <- scas(stks, idxs, fmodel=list(fmod), qmodel=list(qmod), fit="MP", workers=cores)

### output
output.t <- stks + fits
output.t <- FLStocks(c(output.t, FLStocks(om=om.t)))

#------------------------------------------------------------------
# OM = lorenz #------------------------------------------------------------------

## Operating model
om.l <- hke.stk
m(om.l)[] <- m.l
#fit <- sca(om.l, hke.idx, fmodel=fmod, qmodel=qmod, fit="MCMC", mc=mc)
# remove mcscale burnin period
#fit <- burnin(fit, burnin)
fit <- simulate(sca(om.l, hke.idx, fmodel=fmod, qmodel=qmod), it, seed=rs)

om.l <- om.l+fit

## OEM
stk.oem <- om.l
idx.oem <- hke.idx
index(idx.oem) <- index(fit)[[1]]
idx.oem <- FLIndices(idx=idx.oem)

## MP

### current
stk.c <- stk.oem
m(stk.c)[] <- m.c

### chen watanabe
stk.cw <- stk.oem
m(stk.cw)[] <- m.cw

### Then
stk.t <- stk.oem
m(stk.t)[] <- m.t

### Lorenz
stk.l <- stk.oem
m(stk.l)[] <- m.l

### Brodziak
stk.b <- stk.oem
m(stk.b)[] <- m.b

### call sca in parallel
stks <- FLStocks(c=stk.c, cw=stk.cw, t=stk.t, l=stk.l, b=stk.b)
idxs <- list(c=idx.oem, cw=idx.oem, t=idx.oem, l=idx.oem, b=idx.oem)
fits <- scas(stks, idxs, fmodel=list(fmod), qmodel=list(qmod), fit="MP", workers=cores)

### output
output.l <- stks + fits
output.l <- FLStocks(c(output.l, FLStocks(om=om.l)))

#------------------------------------------------------------------
# OM = brodziak #------------------------------------------------------------------

## Operating model
om.b <- hke.stk
m(om.b)[] <- m.b
#fit <- sca(om.b, hke.idx, fmodel=fmod, qmodel=qmod, fit="MCMC", mc=mc)
# remove mcscale burnin period
#fit <- burnin(fit, burnin)
fit <- simulate(sca(om.b, hke.idx, fmodel=fmod, qmodel=qmod), it, seed=rs)

om.b <- om.b+fit

## OEM
stk.oem <- om.b
idx.oem <- hke.idx
index(idx.oem) <- index(fit)[[1]]
idx.oem <- FLIndices(idx=idx.oem)

## MP

### current
stk.c <- stk.oem
m(stk.c)[] <- m.c

### chen watanabe
stk.cw <- stk.oem
m(stk.cw)[] <- m.cw

### Then
stk.t <- stk.oem
m(stk.t)[] <- m.t

### Lorenz
stk.l <- stk.oem
m(stk.l)[] <- m.l

### Brodziak
stk.b <- stk.oem
m(stk.b)[] <- m.b

### call sca in parallel
stks <- FLStocks(c=stk.c, cw=stk.cw, t=stk.t, l=stk.l, b=stk.b)
idxs <- list(c=idx.oem, cw=idx.oem, t=idx.oem, l=idx.oem, b=idx.oem)
fits <- scas(stks, idxs, fmodel=list(fmod), qmodel=list(qmod), fit="MP", workers=cores)

### output
output.b <- stks + fits
output.b <- FLStocks(c(output.b, FLStocks(om=om.b)))

#------------------------------------------------------------------
# Plots of scenarios
#------------------------------------------------------------------

panel_labels <- c("c" = "Current", "cw" = "Chen Whatanabe", "t" = "Then", "l" = "Lorenz", "b" = "Brodziak", "om" = "Operating Model", "Rec" = "Recruitment", "SB" = "Spawning Stock Biomass", "C" = "Catches", "F" = "Fishing mortality")

plot(output.c, probs=c(0.9,0.9,0.5,0.1,0.1)) +
  facet_grid(qname ~ stock, scales = "free_y", labeller = as_labeller(panel_labels)) +
  theme(legend.position = "none")

plot(output.cw, probs=c(0.9,0.9,0.5,0.1,0.1)) +
  facet_grid(qname ~ stock, scales = "free_y", labeller = as_labeller(panel_labels)) +
  theme(legend.position = "none")

plot(output.t, probs=c(0.9,0.9,0.5,0.1,0.1)) +
  facet_grid(qname ~ stock, scales = "free_y", labeller = as_labeller(panel_labels)) +
  theme(legend.position = "none")

plot(output.l, probs=c(0.9,0.9,0.5,0.1,0.1)) +
  facet_grid(qname ~ stock, scales = "free_y", labeller = as_labeller(panel_labels)) +
  theme(legend.position = "none")

plot(output.b, probs=c(0.9,0.9,0.5,0.1,0.1)) +
  facet_grid(qname ~ stock, scales = "free_y", labeller = as_labeller(panel_labels)) +
  theme(legend.position = "none")

#------------------------------------------------------------------
# Check
#------------------------------------------------------------------

flqs0 <- lapply(output.c, "fbar")
df0 <- data.frame(om="c", metric="f", data=c(flqs0$c/flqs0$om))
flqs0 <- lapply(output.c, "ssb")
df0 <- rbind(df0, data.frame(om="c", metric="ssb", data=c(flqs0$c/flqs0$om)))
flqs0 <- lapply(output.c, "rec")
df0 <- rbind(df0, data.frame(om="c", metric="rec", data=c(flqs0$c/flqs0$om)))

flqs0 <- lapply(output.cw, "fbar")
df0 <- rbind(df0, data.frame(om="cw", metric="f", data=c(flqs0$cw/flqs0$om)))
flqs0 <- lapply(output.cw, "ssb")
df0 <- rbind(df0, data.frame(om="cw", metric="ssb", data=c(flqs0$cw/flqs0$om)))
flqs0 <- lapply(output.cw, "rec")
df0 <- rbind(df0, data.frame(om="cw", metric="rec", data=c(flqs0$cw/flqs0$om)))

flqs0 <- lapply(output.t, "fbar")
df0 <- rbind(df0, data.frame(om="t", metric="f", data=c(flqs0$t/flqs0$om)))
flqs0 <- lapply(output.t, "ssb")
df0 <- rbind(df0, data.frame(om="t", metric="ssb", data=c(flqs0$t/flqs0$om)))
flqs0 <- lapply(output.t, "rec")
df0 <- rbind(df0, data.frame(om="t", metric="rec", data=c(flqs0$t/flqs0$om)))

flqs0 <- lapply(output.b, "fbar")
df0 <- rbind(df0, data.frame(om="b", metric="f", data=c(flqs0$b/flqs0$om)))
flqs0 <- lapply(output.b, "ssb")
df0 <- rbind(df0, data.frame(om="b", metric="ssb", data=c(flqs0$b/flqs0$om)))
flqs0 <- lapply(output.b, "rec")
df0 <- rbind(df0, data.frame(om="b", metric="rec", data=c(flqs0$b/flqs0$om)))

flqs0 <- lapply(output.l, "fbar")
df0 <- rbind(df0, data.frame(om="l", metric="f", data=c(flqs0$l/flqs0$om)))
flqs0 <- lapply(output.l, "ssb")
df0 <- rbind(df0, data.frame(om="l", metric="ssb", data=c(flqs0$l/flqs0$om)))
flqs0 <- lapply(output.l, "rec")
df0 <- rbind(df0, data.frame(om="l", metric="rec", data=c(flqs0$l/flqs0$om)))

histogram(~data|om*metric, data=df0)

checkbar.df <- with(df0, tapply(data, list(om,metric), mean))
checksd.df <- with(df0, tapply(data, list(om,metric), sd))
checkmax.df <- with(df0, tapply(data, list(om,metric), max))
checkmin.df <- with(df0, tapply(data, list(om,metric), min))


kable(checkbar.df, digit=3, caption="Mean ratio between EM of OM and OM metrics")

kable(checksd.df, digit=3, caption="Standard deviation of ratio between EM of OM and OM metrics")

kable(checkmax.df, digit=3, caption="Max ratio between EM of OM and OM metrics")

kable(checkmin.df, digit=3, caption="Min ratio between EM of OM and OM metrics")

#------------------------------------------------------------------
# Results and tables
#------------------------------------------------------------------

output.all <- list("Current" = output.c[c("cw", "t","l","b","om")],
  "ChenWatanabe" = output.cw[c("c", "t","l","b","om")],
  "Then"= output.t[c("c", "cw","l","b","om")],
  "Lorenz"= output.l[c("c", "cw", "t","b","om")],
  "Brodziak"= output.b[c("c", "cw", "t","l","om")]
  )

lst0 <- lapply(output.all, function(x){
  # This is not the most efficient code ...
  obj0 <- lapply(x, rec)
  obj0 <- as.data.frame(obj0)
  obj0$om <- obj0[obj0$qname=="om","data"]
  #rec=with(subset(obj0, obj0$qname!="om"), sqrt(mean(data-om)^2))
  rec <- with(subset(obj0, obj0$qname!="om"), mad(data-om))

  obj0 <- lapply(x, ssb)
  obj0 <- as.data.frame(obj0)
  obj0$om <- obj0[obj0$qname=="om","data"]
  #ssb=with(subset(obj0, obj0$qname!="om"), sqrt(mean((data-om)^2))
  ssb <- with(subset(obj0, obj0$qname!="om"), mad(data-om))

  obj0 <- lapply(x, fbar)
  obj0 <- as.data.frame(obj0)
  obj0$om <- obj0[obj0$qname=="om","data"]
  #fbar=with(subset(obj0, obj0$qname!="om"), sqrt(mean((data-om)^2))
  fbar <- with(subset(obj0, obj0$qname!="om"), mad(data-om))

  c(rec, ssb, fbar)
})

results_om <- do.call(cbind, lst0)
results_om[] <- t(apply(results_om, 1, function(x) x/min(x)))
rownames(results_om) <- c("recruitment", "ssb", "f")

kable(results_om, digit=2, caption="MAD for each OM divided my the maximum to simplify the comparison")


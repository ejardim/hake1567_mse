##### HKE IN GSA 1-5-6-7, 2007-2019

rm(list=ls())
graphics.off()

#load the libraries:
library(FLCore)
# library(FLEDA)
# library(FLXSA)
library(FLAssess)
library(FLash)
library(FLa4a)
library(ggplotFL)
library(tidyverse)
library(ggplot2)
library(ggpubr)
require(plyr)
require(FLBRP)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


##########################################################################################
# load the diagnostic functions: retro, mohns_rho, max_retro, min_retro
source("to_be_sourced/diagnostic_functions.R")
# load the estimation functions: selectivity plots (by age and year) and F by age over years
source("to_be_sourced/estimates_functions.R")


######################################################################################
# define the model run number/name and create the output folder
dir.create(file.path("../output/"))

out_folder <- file.path("../output/")
out_indat_fold <- file.path(out_folder, "plots_input_data")
dir.create(out_indat_fold)

out_diag_fold <- file.path(out_folder, "diagnostics")
dir.create(out_diag_fold)

out_est_fold <- file.path(out_folder, "model_estimates")
dir.create(out_est_fold)

out_obj_fold <- file.path(out_folder, "fit_model_objects")
dir.create(out_obj_fold)



###################################################################
#read stock and idx files
hke.stk <- readRDS(file="../../data_preparation/output/final_idx_stk_objects/HKE_1_5_6_7_stk.rds")
hke.idx <- readRDS(file="../../data_preparation/output/final_idx_stk_objects/HKE_1_5_6_7_idx.rds")


##############################################################################################################
## Plot input data

#pdf("./diagnostics_hake1567.pdf")
jpeg(file.path(out_indat_fold, "catch_at_age.jpg"), res=550, height=3000, width=4500)
bubbles(age~year, data=(catch.n(hke.stk)), bub.scale=10)
dev.off()

jpeg(file.path(out_indat_fold, "Stock_weigth_at_age.jpg"), res=550, height=3000, width=4500)
bubbles(age~year, data=(catch.wt(hke.stk)), bub.scale=5)
dev.off()

jpeg(file.path(out_indat_fold, "IDX_age_distr.jpg"), res=550, height=3000, width=4500)
bubbles(age~year, data=(index(hke.idx)), bub.scale=5)   
dev.off()

#plot the stock
jpeg(file.path(out_indat_fold, "Observed_catch.jpg"), res=550, height=3000, width=4500)
plot(hke.stk)   ## 
dev.off()

jpeg(file.path(out_indat_fold, "Observed_catch_perAge.jpg"), res=550, height=3000, width=4500)
plot(catch.n(hke.stk))
dev.off()

jpeg(file.path(out_indat_fold, "Observed_survey_idx_perAge.jpg"), res=550, height=3000, width=4500)
plot(index(hke.idx))
dev.off()

jpeg(file.path(out_indat_fold, "weight_at_age.jpg"), res=550, height=3000, width=4500)
plot(catch.wt(hke.stk))
dev.off()


## plotting cohorts
#call the function to produce the plots of cohort consistency
source("to_be_sourced/Cohort_Consistency.r")
#cccor(pil.tun0) ## three options: "pearson","kendall" or "spearman"
cccor(hke.stk,"pearson","HKE_1567")
cccor(FLIndices(hke.idx),"pearson", "HKE_1567")

# To plot numbers at age by year
jpeg(file.path(out_indat_fold, "C_at_age.jpg"), res=450, height=3000, width=4500)
xyplot(data~age, group=year, data=catch.n(hke.stk), main="Catches age structure", type=c("g", "l"), ylab= "N (thousands)", auto.key=list(title="Year",points=F, lines=T, space="right"))
dev.off()

jpeg(file.path(out_indat_fold, "Index_at_age.jpg"), res=450, height=3000, width=4500)
xyplot(data~age, group=year, data=index(hke.idx), main="Survey age structure", type=c("g", "l"), ylab= "", auto.key=list(title="Year",points=F, lines=T, space="right"))
dev.off()




###################################################################
############## DIFFERENT QMODEL, FMODEL and SRMODEL

# Survey catchability submodel
qmodel<- list(~ I(1/(1+exp(-age)))) ###

# Stock-Recruitment submodels
srmodel <- ~ factor(year)

# Fishing mortality submodels
fmodel157 <- ~s(age, k = 4) + s(year, k = 8) + te(age, year, k = c(3, 10))


#====================================================================
# Fit the assessment models
#====================================================================

fit <- sca(hke.stk, hke.idx, fmodel=fmodel157, qmodel=qmodel, srmodel=srmodel)

# save(fit, file=paste0("../objects/", model_run, "_HKE_1_5_6_7_fit.RData"))
saveRDS(fit, "../output/fit_model_objects/HKE_1_5_6_7_FIT.rds")


stk <- hke.stk + fit

# save(stk, file=paste0("../objects/", model_run, "_HKE_1_5_6_7_model.RData"))
saveRDS(stk, "../output/fit_model_objects/HKE_1_5_6_7_MODEL.rds")

###################
# model estimates
png(file.path(out_est_fold, "model_estimates.png"), height=2500, width=3500, res=350, pointsize = 12)
print(plot(stk))
dev.off()

# plot estimated Q at age and year
# compute N for the fraction of the year the survey is carried out
sfrac <- mean(range(hke.idx)[c("startf", "endf")])
# fraction of total mortality up to that moment
Z <- (m(hke.stk) + harvest(fit))*sfrac
lst <- dimnames(fit@index)
# survivors
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn[-6,,,,,,]
png(file.path(out_est_fold, "fitted_Q-at-age-and-year.png"), height=2500, width=3500, res=350, pointsize = 12)
print(wireframe(qhat, zlab="Q",screen = list(z = 30, x = -60)))
dev.off()

png(file.path(out_est_fold, "/F_wireframe.png"), res=450, height=3000, width=4500)
wireframe(data ~ age + year, data = as.data.frame(harvest(fit)), zlab="F", drape = TRUE, screen = list(x = -40, y=-65, z=-50))
dev.off()


# Selectivity
png(file.path(out_est_fold, "Selectivity_by_age.png"), height=2500, width=3500, res=350, pointsize = 12)
print(plotSels_by_age(hke.stk + fit))
dev.off()

png(file.path(out_est_fold, "Selectivity_by_year.png"), height=2500, width=3500, res=350, pointsize = 12)
print(plotSels_by_year(hke.stk + fit))
dev.off()

# F curves
png(file.path(out_est_fold, "Fs_by_age.png"), height=2500, width=3500, res=350, pointsize = 12)
print(plotFs_by_age(hke.stk + fit))
dev.off()

# Save some model estimates
Fbar <- as.data.frame(fbar(stk))
write.csv(Fbar, file.path(out_est_fold, "Fbar.csv"))

SSB <- as.data.frame(ssb(stk))
write.csv(SSB, file.path(out_est_fold, "SSB.csv"))





##########################
# Diagnostics

res <- residuals(fit, hke.stk, hke.idx)

# fit summary
fitsum <- as.data.frame(fitSumm(fit))
colnames(fitsum) <- c("value")
fitsum <- fitsum %>% 
  mutate(parameter=rownames(.)) %>% 
  relocate(parameter, value)
write.csv(fitsum, file=file.path(out_diag_fold, "fit_summary.csv"), row.names=FALSE)


png(file.path(out_diag_fold, "residuals_totalcatch.png"), height=2500, width=3500, res=350, pointsize = 12)
tot_catch_res<-computeCatchDiagnostics(fit, hke.stk)
tot_catch_res_plot<-plot(tot_catch_res)
plot(tot_catch_res_plot[3])
dev.off()

png(file.path(out_diag_fold, "residuals_lineplot.png"), height=2500, width=3500, res=350, pointsize = 12)
print(plot(res))
dev.off()

png(file.path(out_diag_fold, "residuals_bubbleplot.png"), height=2500, width=3500, res=350, pointsize = 12)
print(bubbles(res))
dev.off()

png(file.path(out_diag_fold, "residuals_qqplot.png"), height=2500, width=3500, res=350, pointsize = 12)
print(qqmath(res))
dev.off()


idx_obs <- FLQuants(index(hke.idx))
idx_obs@names <- "IND"

png(file.path(out_diag_fold, "combined_residuals.png"), height=2500, width=3500, res=350, pointsize = 12)
plot_index_combined<-plotRunstest(index(fit) , idx_obs, combine=T)+
  theme_bw() +ylab("Index")+ggtitle("Combined Residulas")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.margin = margin(5, 5, 0, 5),
        plot.title = element_text(hjust = 0.5),strip.text = element_blank(),panel.spacing = unit(0, "lines"),
        legend.position = "none")

plot_catch_combined<-plotRunstest(catch.n(fit),catch.n(hke.stk), combine = T)+
  theme_bw() +ylab("Catch")+
  theme(plot.margin = margin(0, 5, 5, 5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_blank(),panel.spacing = unit(0, "lines"),legend.position = "none")
ggarrange(plot_index_combined,plot_catch_combined, nrow=2, align = "v")
dev.off()



png(file.path(out_diag_fold, "combined_residuals_by_age.png"), height=2500, width=3500, res=350, pointsize = 12)

plot_age_index<-plotRunstest(index(fit) , idx_obs, combine=F)+
  theme_bw() + facet_wrap(~age,nrow=1)+ylab("Index")+ggtitle("Residulas by age")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.margin = margin(5, 5, 0, 5),
        plot.title = element_text(hjust = 0.5),panel.spacing = unit(0, "lines"),legend.position = "none")

plot_age_catch<-plotRunstest(catch.n(fit),catch.n(hke.stk), combine = F)+
  theme_bw() + facet_wrap(~age,nrow=1)+ylab("Catch")+
  theme(plot.margin = margin(0, 5, 5, 5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing = unit(0, "lines"),legend.position = "none")
ggarrange(plot_age_index,plot_age_catch, nrow=2, align = "v")

dev.off()


# Fit to observed data
png(file.path(out_diag_fold, "fitted_catch-at-age.png"), height=2500, width=3500, res=350, pointsize = 12)
print(plot(fit, hke.stk))
dev.off()

png(file.path(out_diag_fold, "fitted_index-at-age.png"), height=2500, width=3500, res=350, pointsize = 12)
print(plot(fit, FLIndices(hke.idx)))
dev.off()

# plot estimated F at age and year
png(file.path(out_diag_fold, "fitted_F-at-age-and-year.png"), height=2500, width=3500, res=350, pointsize = 12)
print(wireframe(harvest(fit), zlab="F",screen = list(z = 30, x = -30)))
dev.off()


# plot est/obs catch at age over time
ca_fit <- as.data.frame(fit@catch.n)
ca_fit$origin <- "fit"
ca_obs <- as.data.frame(hke.stk@catch.n)
ca_obs$origin <- "obs"
ca <- rbind(ca_fit, ca_obs)

png(filename=file.path(out_diag_fold, "obs_est_catch_age_year.png"), height=3500, width=5500, res=500)
ggplot2::ggplot() +
  geom_line(data=ca_fit, aes(x=year, y=data), col="blue") +
  geom_point(data=ca_obs, aes(x=year, y=data), col="black") +
  facet_wrap(~age, scales="free")
dev.off()

# plot est/obs survey index at age over time
idx_fit <- as.data.frame(fit@index)
idx_fit$origin <- "fit"
idx_fit <- idx_fit[,!colnames(idx_fit) %in% c("qname")]
idx_obs <- as.data.frame(hke.idx@index)
idx_obs$origin <- "obs"
idx <- rbind(idx_fit, idx_obs)

png(filename=file.path(out_diag_fold, "obs_est_survey_idx_age_year.png"), height=3500, width=5500, res=500)
ggplot2::ggplot() +
  geom_line(data=idx_fit, aes(x=year, y=data), col="blue") +
  geom_point(data=idx_obs, aes(x=year, y=data), col="black") +
  facet_wrap(~age, scales="free")
dev.off()



# #####################
# # MASE
# # make xval object to run and plot retro and hindcast
# # set nyears for retrospective (the choice should be done base on the available number of years along the time series)
# nyears <- 5
# # set number of year for average biology and selectivity
# nsq <- 3
# 
# hc <- a4ahcxval(stk, FLIndices(hke.idx),
#                 nyears = nyears, nsq = nsq, srmodel = srmodel(fit), fmodel = fmodel(fit),
#                 qmodel = formula(qmodel(fit)))
# 
# MASE <- mase(idx, hc$indices)
# 
# # save MASE results for diags table
# MASE_result <- data.frame(test = c("mase"), value = as.numeric(c(MASE)))%>%
#   mutate(interpretation = case_when(value < 1 ~ "good", value > 1 ~ "check"))




# Retrospective pattern
# hc <- a4ahcxval(hke.stk, hke.idx, nyears = nyears, nsq = nsq, srmodel = srmodel(fit), fmodel = fmodel(fit), qmodel = formula(qmodel(fit)))
# plot(hc$stocks)
# mohn_ssb <- icesAdvice::mohn(mohnMatrix(hc$stocks, ssb))
# mohn_f <- icesAdvice::mohn(mohnMatrix(hc$stocks, fbar))
# 

# Retrospective pattern

# Create the retro function
retro <- function(stk, idxs, retro=5, kfrac="missing", k, k2, ...){
  args <- list(...)
  if(missing(kfrac)) kfrac <- 0.3
  lst0 <- split(0:retro, 0:retro)
  lst0 <- lapply(lst0, function(x){
    yr <- range(stk)["maxyear"] - x
    args$stock <- window(stk, end=yr)
    args$indices <- FLIndices(window(idxs, end=yr))
    KY <- unname(k - floor(x*kfrac))
    KZ <- unname(k2 - floor(x*kfrac))
    # NOTE THIS HAS TO BE ADAPTED FOR THE MODEL USED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fmod <- substitute(~s(age, k = 4) + s(year, k = KY) + te(age, year, k = c(3, 10)), list(KY=KY))
    args$fmodel <- as.formula(fmod)
    # args$srmodel <- as.formula(srmod)
    fit <- do.call("sca", args)
    args$stock + fit
  })
  FLStocks(lst0)
}
retro.lst <- retro(hke.stk, hke.idx, kfrac=0.3, k=8, k2=NA, retro=5, qmodel=qmodel) # The only submodel that needs to be indicated is the qmodel. The fmodel and srmodel are taken from the fit object.
mohn_rho <- mohn(retro.lst) 
png(file.path(out_diag_fold, "Retrospective.png"), res=450, height=3000, width=4500)
print(plot(retro.lst))
dev.off()
write.csv(mohn_rho, file.path(out_diag_fold, "mohn_rho.csv"), row.names = FALSE)

# max and min retro
max_retro <- maxretro(retro.lst)
write.csv(max_retro, file.path(out_diag_fold, "max_retro.csv"), row.names = FALSE)
min_retro <- minretro(retro.lst)
write.csv(min_retro, file.path(out_diag_fold, "min_retro.csv"), row.names = FALSE)




#######################################
########### Simulate 

# estimate confidence intervals for model estimates by using de var-covar matrix from hessian

fit.sim <- simulate(fit, 1000)

# update stock object
stk.sim <- hke.stk+fit.sim

# Retro pattern with confidence intervals
png(file.path(out_diag_fold, "Retrospective_with_conf_interval_stksim.png"), res=450, height=3000, width=4500)
print(plot(stk.sim, retro.lst))
dev.off()


#iterVars(catch.n(fit.sim))

jpeg(file.path(out_est_fold, "simulated_summary.jpg"), res=550, height=3000, width=4500)
plot(stk.sim)
dev.off()

jpeg(file.path(out_est_fold, "fbar.jpg"), res=550, height=3000, width=4500)
plot(fbar(stk.sim))
dev.off()

jpeg(file.path(out_est_fold, "ssb.jpg"), res=550, height=3000, width=4500)
plot(ssb(stk.sim))
dev.off()

jpeg(file.path(out_est_fold, "recr.jpg"), res=550, height=3000, width=4500)
plot(rec(stk.sim))
dev.off()

jpeg(file.path(out_est_fold, "catch.jpg"), res=550, height=3000, width=4500)
plot(catch(stk.sim))
dev.off()


### REtro with confidence interval from simulations. Done in a rush, this part should be better coded to reduce extension
library(dplyr)
tmp <- as.data.frame(ssb(stk.sim)) %>%
  group_by(year) %>%
  dplyr::summarize(mean=mean(data, na.rm=TRUE),
            phigh=quantile(data, prob=0.975, na.rm=TRUE),
            plow=quantile(data, prob=0.025, na.rm=TRUE))
tmp$var <- "ssb"

ssb_retro <- matrix(NA, nrow=length(c(dims(stk.sim)$minyear:dims(stk.sim)$maxyear)), ncol=6)
ssb_retro[,1] <- ssb(retro.lst[[1]])
ssb_retro[,2] <- c(ssb(retro.lst[[2]]), NA)
ssb_retro[,3] <- c(ssb(retro.lst[[3]]), NA, NA)
ssb_retro[,4] <- c(ssb(retro.lst[[4]]), NA, NA, NA)
ssb_retro[,5] <- c(ssb(retro.lst[[5]]), NA, NA, NA, NA)
ssb_retro[,6] <- c(ssb(retro.lst[[6]]), NA, NA, NA, NA, NA)

  ssb_retro <- data.frame(ssb_retro)
colnames(ssb_retro) <- c("retro_0", "retro_1", "retro_2", "retro_3", "retro_4", "retro_5")

datssb <- cbind(tmp, ssb_retro)




tmp <- as.data.frame(catch(stk.sim)) %>%
  group_by(year) %>%
  dplyr::summarize(mean=mean(data, na.rm=TRUE),
            phigh=quantile(data, prob=0.975, na.rm=TRUE),
            plow=quantile(data, prob=0.025, na.rm=TRUE))
tmp$var <- "catch"

catch_retro <- matrix(NA, nrow=length(c(dims(stk.sim)$minyear:dims(stk.sim)$maxyear)), ncol=6)
catch_retro[,1] <- catch(retro.lst[[1]])
catch_retro[,2] <- c(catch(retro.lst[[2]]), NA)
catch_retro[,3] <- c(catch(retro.lst[[3]]), NA, NA)
catch_retro[,4] <- c(catch(retro.lst[[4]]), NA, NA, NA)
catch_retro[,5] <- c(catch(retro.lst[[5]]), NA, NA, NA, NA)
catch_retro[,6] <- c(catch(retro.lst[[6]]), NA, NA, NA, NA, NA)

catch_retro <- data.frame(catch_retro)
colnames(catch_retro) <- c("retro_0", "retro_1", "retro_2", "retro_3", "retro_4", "retro_5")

datcatch <- cbind(tmp, catch_retro)


tmp <- as.data.frame(fbar(stk.sim)) %>%
  group_by(year) %>%
  dplyr::summarize(mean=mean(data, na.rm=TRUE),
            phigh=quantile(data, prob=0.975, na.rm=TRUE),
            plow=quantile(data, prob=0.025, na.rm=TRUE))
tmp$var <- "fbar"

fbar_retro <- matrix(NA, nrow=length(c(dims(stk.sim)$minyear:dims(stk.sim)$maxyear)), ncol=6)
fbar_retro[,1] <- fbar(retro.lst[[1]])
fbar_retro[,2] <- c(fbar(retro.lst[[2]]), NA)
fbar_retro[,3] <- c(fbar(retro.lst[[3]]), NA, NA)
fbar_retro[,4] <- c(fbar(retro.lst[[4]]), NA, NA, NA)
fbar_retro[,5] <- c(fbar(retro.lst[[5]]), NA, NA, NA, NA)
fbar_retro[,6] <- c(fbar(retro.lst[[6]]), NA, NA, NA, NA, NA)

fbar_retro <- data.frame(fbar_retro)
colnames(fbar_retro) <- c("retro_0", "retro_1", "retro_2", "retro_3", "retro_4", "retro_5")

datfbar <- cbind(tmp, fbar_retro)




tmp <- as.data.frame(rec(stk.sim)) %>%
  group_by(year) %>%
  dplyr::summarize(mean=mean(data, na.rm=TRUE),
            phigh=quantile(data, prob=0.975, na.rm=TRUE),
            plow=quantile(data, prob=0.025, na.rm=TRUE))
tmp$var <- "rec"

rec_retro <- matrix(NA, nrow=length(c(dims(stk.sim)$minyear:dims(stk.sim)$maxyear)), ncol=6)
rec_retro[,1] <- rec(retro.lst[[1]])
rec_retro[,2] <- c(rec(retro.lst[[2]]), NA)
rec_retro[,3] <- c(rec(retro.lst[[3]]), NA, NA)
rec_retro[,4] <- c(rec(retro.lst[[4]]), NA, NA, NA)
rec_retro[,5] <- c(rec(retro.lst[[5]]), NA, NA, NA, NA)
rec_retro[,6] <- c(rec(retro.lst[[6]]), NA, NA, NA, NA, NA)

rec_retro <- data.frame(rec_retro)
colnames(rec_retro) <- c("retro_0", "retro_1", "retro_2", "retro_3", "retro_4", "retro_5")

datrec <- cbind(tmp, rec_retro)

dat <- rbind(datssb, datcatch, datfbar, datrec)


p <- dat %>% 
  ggplot(., aes(x=year, y=retro_0)) + 
  geom_line() +
  geom_ribbon(aes(ymin=plow, ymax=phigh, alpha=0.5)) +
  geom_line(aes(y=retro_1, col="retro_1")) +
  geom_line(aes(y=retro_2, col="retro_2")) +
  geom_line(aes(y=retro_3, col="retro_3")) +
  geom_line(aes(y=retro_4, col="retro_4")) +
  geom_line(aes(y=retro_5, col="retro_5")) +
  facet_wrap(~var, scales="free")

png(filename=file.path(out_diag_fold, "retro_with_conf_int.png"), height=3000, width=4500, res=500, pointsize = 12)
print(p)            
dev.off()



#### Plot estimations

# F at age
harvest(stk)
fatage <- as.data.frame(harvest(stk)) %>% 
  dplyr::select(year, age, data) %>% 
  tidyr::pivot_wider(names_from = year, values_from = data)
write.csv(fatage, file.path(out_est_fold, "f_at_age_year.csv"), row.names = F)

png(filename=file.path(out_est_fold, "harvest_at_age.png"), height=3000, width=4500, res=500, pointsize = 12)
plot(harvest(stk))
dev.off()


png(filename=file.path(out_est_fold, "stock_n_at_age.png"), height=3000, width=4500, res=500, pointsize = 12)
plot(stk@stock.n)
dev.off()

stock.n(stk)
natage <- as.data.frame(stock.n(stk)) %>% 
  dplyr::select(year, age, data) %>% 
  tidyr::pivot_wider(names_from = year, values_from = data)
write.csv(natage, file.path(out_est_fold, "n_at_age_year.csv"), row.names = F)


png(filename=file.path(out_est_fold, "catch_at_age.png"), height=3000, width=4500, res=500, pointsize = 12)
plot(stk@catch.n)
dev.off()

catch.n(stk)
catage <- as.data.frame(catch.n(stk)) %>% 
  dplyr::select(year, age, data) %>% 
  tidyr::pivot_wider(names_from = year, values_from = data)
write.csv(catage, file.path(out_est_fold, "n_at_age_year.csv"), row.names = F)



# #stk <- stk_rec
# stk_brp <- brp(FLBRP(stk))
# refpts(stk_brp)
# f01 <- c(refpts(stk_brp)["f0.1","harvest"])
# f01
# #0.4058128
# fbar(stk)
# #1.3176
# 
# fbar(stk)[,"2022"]/f01
# #3.2469
# #
# 
# 
# plot(ypr(stk_brp)~fbar(stk_brp),type='l')
# 
# # jpeg("../output/output_assess/Repts_plot_old.jpg", res=450, height=3000, width=4500)
# plot(stk_brp)
# # dev.off()


otp <- cbind(round(an(catch(stk)),0), round(an(ssb(stk)),0), round(an(rec(stk)),0), round(an(fbar(stk)),3))
colnames(otp) <-c("catch", "ssb", "recruitment", "Fbar")
otp <- data.frame(otp)
dir.create("../output/Report/")
write.csv(otp, file="../output/Report/stk_output.csv",row.names = FALSE)
write.csv(otp, file=file.path(out_est_fold, "stk_output.csv"),row.names = FALSE)


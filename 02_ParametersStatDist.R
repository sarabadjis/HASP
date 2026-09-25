## Processing ATUS Data
# August 2026
# Sophia Arabadjis

library(here)
library(tidyverse)
library(survival)
library(survminer)
library(flexsurv)
library(gamlss.dist)
library(fitdistrplus)
# Call instead of load
#library(fitdistrplus)
#library(hms)

source(here("code","00_Source.R"))

# ---- Read in Data ----

diary_ftc <- readRDS(here("data/proc_data","01_FTS.rds"))
diary_ptc <- readRDS(here("data/proc_data","01_PTS.rds"))

## ---- Quick Visuals ----

#Time to meals
surv_plots.f("Eating and Drinking",n=1)
surv_plots.f("Eating and Drinking",n=2)
surv_plots.f("Eating and Drinking",n=3)

# Time to sleeping episodes
surv_plots.f("Sleeping",n=1)
surv_plots.f("Sleeping",n=2)

## ---- Fitting Sleep Distributions ----

# Sleep Densities
sleepdensftc <- diary_ftc %>%
  group_by(TUCASEID,activity)%>%
  summarise(mins = sum(TUACTDUR24))%>%
  ungroup()%>%
  dplyr::filter(activity == "Sleeping")%>%
  dplyr::pull(mins)%>%
  density(window = "e")

sleepdensptc <- diary_ptc %>%
  group_by(TUCASEID,activity)%>%
  summarise(mins = sum(TUACTDUR24))%>%
  ungroup()%>%
  dplyr::filter(activity == "Sleeping")%>%
  dplyr::pull(mins)%>%
  density(window = "e")

# Sleep Data
sleepdatftc <- diary_ftc %>%
  group_by(TUCASEID,activity)%>%
  summarise(mins = sum(TUACTDUR24))%>%
  ungroup()%>%
  dplyr::filter(activity == "Sleeping")%>%
  dplyr::pull(mins)
sleepdatptc <- diary_ptc %>%
  group_by(TUCASEID,activity)%>%
  summarise(mins = sum(TUACTDUR24))%>%
  ungroup()%>%
  dplyr::filter(activity == "Sleeping")%>%
  dplyr::pull(mins)

hist(sleepdatftc,breaks = 15,main = "Histogram Sleep Durations\n Full-Time Students",
     xlab = "Minutes")
hist(sleepdatptc,breaks = 15, main = "Histogram Sleep Durations\n Part-Time Students",
     xlab = "Minutes")

fit_n_ftc <- fitdistrplus::fitdist(sleepdatftc,distr = "norm")
fit_w_ftc <- fitdistrplus::fitdist(sleepdatftc,distr = "weibull")

fit_n_ptc <- fitdistrplus::fitdist(sleepdatptc,distr = "norm")
fit_w_ptc <- fitdistrplus::fitdist(sleepdatptc,distr = "weibull")
fit_g_ptc <- fitdistrplus::fitdist(sleepdatptc,distr = "gamma")

fitdistrplus::cdfcomp(list(fit_n_ptc,fit_w_ptc,fit_g_ptc))
# Per AIC, weibull fits best

params_ftc <- coef(fit_n_ftc)
params_ptc <- coef(fit_w_ptc)

# Sim Check

s_ptc <- replicate(100,rweibull(500,shape = params_ptc[[1]],scale = params_ptc[[2]]))
d_ptc <- apply(X = s_ptc,MARGIN = 2,FUN = density,window="e")

par(mfrow = c(1, 2),oma = c(0, 1, 3, 0))
plot(sleepdensptc,ylim = c(0,.005), 
     main = "Part-Time Students\n(Weibull)",
     bty='n',ylab = "",las=1,cex=.8)
for(i in 1:100){
  lines(d_ptc[[i]],col='grey')
}
lines(sleepdensptc)

s_ftc <- replicate(100,rnorm(500,mean = params_ftc[[1]],sd = params_ftc[[2]]))
d_ftc <- apply(X = s_ftc,MARGIN = 2,FUN = density, window="e")

plot(sleepdensftc,ylim = c(0,.005), 
     main = "Full-Time Students\n(Gaussian)",
     bty='n',las=1, ylab = "")
for(i in 1:100){
  lines(d_ftc[[i]],col='grey')
}
lines(sleepdensftc)
mtext("Simulated Density Check: Sleep Durations", side = 3, outer = TRUE, cex = 1.5, font = 2)
mtext("Density", side = 2, outer = TRUE, cex = 1, line = -.1, font = 1)
dev.off()


## ---- Distributional Fits for Meal Timing ----

m_ft <- diary_ftc %>%
  group_by(TUCASEID,activity)%>%
  arrange(as.numeric(TUACTIVITY_N))%>%
  mutate(mnum = ifelse(activity == "Eating and Drinking",row_number(),0),
         snum = ifelse(activity == "Sleeping",row_number(),0))%>%
  ungroup()

m_pt <- diary_ptc %>%
  group_by(TUCASEID,activity)%>%
  arrange(as.numeric(TUACTIVITY_N))%>%
  mutate(mnum = ifelse(activity == "Eating and Drinking",row_number(),0),
         snum = ifelse(activity == "Sleeping",row_number(),0))%>%
  ungroup()

# Time to first meal
m1 <- m_ft %>% 
  mutate(activity = case_when(mnum >= 1 ~ paste("Meal",mnum,sep="_"),
                              snum >= 1 ~ paste("Sleeping",snum,sep="_"), 
                              TRUE ~ activity)) %>%
  dplyr::filter(activity %in% c("Meal_1","Sleeping_1"))%>%
  dplyr::select(TUCASEID,activity,startact,stopact)%>%
  pivot_wider(id_cols = TUCASEID,names_from = activity,values_from = c(startact,stopact))%>%
  mutate(meal1dist = hms::as_hms(as.numeric(startact_Meal_1-stopact_Sleeping_1)),
         meal1dist = ifelse(sign(meal1dist)== -1,meal1dist+hms::hms(seconds = 0, minutes = 0, hours=24),meal1dist),
         meal1dist = hms::hms(seconds = meal1dist))%>%
  mutate(m1d = as.numeric(meal1dist)/60)%>%
  dplyr::filter(!is.na(m1d))

m1pt <- m_pt %>% 
  mutate(activity = case_when(mnum >= 1 ~ paste("Meal",mnum,sep="_"),
                              snum >= 1 ~ paste("Sleeping",snum,sep="_"), 
                              TRUE ~ activity)) %>%
  dplyr::filter(activity %in% c("Meal_1","Sleeping_1"))%>%
  dplyr::select(TUCASEID,activity,startact,stopact)%>%
  pivot_wider(id_cols = TUCASEID,names_from = activity,values_from = c(startact,stopact))%>%
  mutate(meal1dist = hms::as_hms(as.numeric(startact_Meal_1-stopact_Sleeping_1)),
         meal1dist = ifelse(sign(meal1dist)== -1,meal1dist+hms::hms(seconds = 0, minutes = 0, hours=24),meal1dist),
         meal1dist = hms::hms(seconds = meal1dist))%>%
  mutate(m1d = as.numeric(meal1dist)/60)%>%
  dplyr::filter(!is.na(m1d))

m1 %>%
  mutate(immediate = ifelse(m1d==0,"Immediate Meal","Delayed Meal"))%>%
  group_by(immediate)%>%
  summarize(n=n())
#32/207 full-time students have immediate meal

m1pt %>%
  mutate(immediate = ifelse(m1d==0,"Immediate Meal","Delayed Meal"))%>%
  group_by(immediate)%>%
  summarize(n=n())
#6/53 part-time students have immediate meal

# ---- Meal 1 Model Specification Search ----

meal1_ft1 <- fitdistrplus::fitdist(data = m1$m1d,distr = "exp")
meal1_ft2 <- fitdistrplus::fitdist(data = m1$m1d,distr = "norm")
meal1_ft3 <- fitdistrplus::fitdist(m1$m1d,distr = "gompertz", start = list(shape = 0.0001,rate = 0.0001))
meal1_ft4 <- fitdistrplus::fitdist(m1$m1d,distr = "ZIP", start = list(mu = 175, sigma = 0.14), discrete = T,
                                   lower=c(0,0), upper = c(Inf, 1))
meal1_ft5 <- fitdistrplus::fitdist(m1$m1d,distr = "ZAZIPF", start = list(mu = 175, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1))
meal1_ft6 <- fitdistrplus::fitdist(m1$m1d,distr = "ZANBI", start = list(mu = 175, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1))
meal1_ft7 <- fitdistrplus::fitdist(m1$m1d,distr = "ZABNB", start = list(mu = 175, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1),method = "mle")
meal1_ft8 <- fitdistrplus::fitdist(m1$m1d,distr = "ZINBI", start = list(mu = 175, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1))
meal1_ft9 <- fitdistrplus::fitdist(m1$m1d,distr = "ZIBNB", start = list(mu = 175, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1),method = "mle")
meal1_ft10 <- fitdistrplus::fitdist(m1$m1d,distr = dZIE, start = list(rate = 0.3,pi = 0.15), discrete=F,
                                   lower=c(rate = 0, pi = 0), upper = c(rate=Inf,pi= 1),method = "mle")
meal1_ft11 <- fitdistrplus::fitdist(m1$m1d,distr = dZIN, start = list(mean = 175, sd=150, pi = 0.1), discrete=F,
                                    lower=c(mean = 0,sd = 0.001,pi = 0), upper = c(mean = Inf, sd = Inf, pi = 1),method = "mle")


fitdistrplus::cdfcomp(list(meal1_ft1,meal1_ft2,meal1_ft3,meal1_ft4, meal1_ft5,meal1_ft6,meal1_ft7,meal1_ft8,meal1_ft9,meal1_ft10,meal1_ft11))
# ZIE is best... 

meal1_pt1 <- fitdistrplus::fitdist(data = m1pt$m1d,distr = "exp")
meal1_pt3 <- fitdistrplus::fitdist(m1pt$m1d,distr = "gompertz", start = list(shape = 0.0001,rate = 0.0001))
meal1_pt4 <- fitdistrplus::fitdist(m1pt$m1d,distr = "ZIP", start = list(mu = 168, sigma = 0.11), discrete = T,
                                   lower=c(0,0), upper = c(Inf, 1))
meal1_pt5 <- fitdistrplus::fitdist(m1pt$m1d,distr = "ZAZIPF", start = list(mu = 168, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1))
meal1_pt6 <- fitdistrplus::fitdist(m1pt$m1d,distr = "ZANBI", start = list(mu = 168, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1))
meal1_pt7 <- fitdistrplus::fitdist(m1pt$m1d,distr = "ZABNB", start = list(mu = 168, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1),method = "mle")
meal1_pt8 <- fitdistrplus::fitdist(m1pt$m1d,distr = "ZINBI", start = list(mu = 168, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1))
meal1_pt9 <- fitdistrplus::fitdist(m1pt$m1d,distr = "ZIBNB", start = list(mu = 168, sigma = .1), discrete=T,
                                   lower=c(0,0), upper = c(Inf, 1),method = "mle")
meal1_pt10 <- fitdistrplus::fitdist(m1pt$m1d,distr = dZIE, start = list(rate = 0.3,pi = 0.11), discrete=T,
                                    lower=c(0,0), upper = c(Inf, 1),method = "mle")
fitdistrplus::cdfcomp(list(meal1_pt1, meal1_pt3, meal1_pt4, meal1_pt5,meal1_pt6,meal1_pt7,meal1_pt8,meal1_pt9,meal1_pt10))

# Extract Coefficient 
Mparams_ptc <- coef(meal1_pt10)
Mparams_ftc <- coef(meal1_ft10)

# Sim Check

s_ptc <- replicate(100,rZIE(500,rate = Mparams_ptc[[1]], pi = Mparams_ptc[[2]]))
d_ptc <- apply(X = s_ptc,MARGIN = 2,FUN = density,window="e")

meal1denspt <- density(m1pt$m1d)

par(mfrow = c(1, 2),oma = c(0, 1, 3, 0))
plot(meal1denspt,ylim = c(0,.005), 
     main = "Part-Time Students\n(Zero-Inflated Exponential)",
     bty='n',ylab = "",las=1,cex=.8)
for(i in 1:100){
  lines(d_ptc[[i]],col='grey')
}
lines(meal1denspt)

meal1densft <- density(m1$m1d)

s_ftc <- replicate(100,rZIE(500,rate = Mparams_ftc[[1]],pi = Mparams_ftc[[2]]))
d_ftc <- apply(X = s_ftc,MARGIN = 2,FUN = density, window="e")

plot(meal1densft,ylim = c(0,.005), 
     main = "Full-Time Students\n(Zero-Inflated Exponential)",
     bty='n',las=1, ylab = "")
for(i in 1:100){
  lines(d_ftc[[i]],col='grey')
}
lines(meal1densft)
mtext("Simulated Density Check: Time-to-First-Meal-From-Waking", side = 3, outer = TRUE, cex = 1.5, font = 2)
mtext("Density", side = 2, outer = TRUE, cex = 1, line = -.1, font = 1)
dev.off()

# ---- Meal 2 Model Specification Search ----

m2 <- m_ft %>% 
  mutate(activity = case_when(mnum >= 1 ~ paste("Meal",mnum,sep="_"),
                              snum >= 1 ~ paste("Sleeping",snum,sep="_"), 
                              TRUE ~ activity)) %>%
  dplyr::filter(activity %in% c("Meal_1","Meal_2")) %>%
  dplyr::select(TUCASEID,activity,startact,stopact)%>%
  pivot_wider(id_cols = TUCASEID,names_from = activity,values_from = c(startact,stopact))%>%
  mutate(meal2dist = hms::as_hms(as.numeric(startact_Meal_2-stopact_Meal_1)),
         meal2dist = ifelse(sign(meal2dist)== -1,meal2dist+hms::hms(seconds = 0, minutes = 0, hours=24),meal2dist),
         meal2dist = hms::hms(seconds = meal2dist))%>%
  mutate(m2d = as.numeric(meal2dist)/60)%>%
  dplyr::filter(!is.na(m2d))

m2pt <- m_pt %>% 
  mutate(activity = case_when(mnum >= 1 ~ paste("Meal",mnum,sep="_"),
                              snum >= 1 ~ paste("Sleeping",snum,sep="_"), 
                              TRUE ~ activity)) %>%
  dplyr::filter(activity %in% c("Meal_1","Meal_2")) %>%
  dplyr::select(TUCASEID,activity,startact,stopact)%>%
  pivot_wider(id_cols = TUCASEID,names_from = activity,values_from = c(startact,stopact))%>%
  mutate(meal2dist = hms::as_hms(as.numeric(startact_Meal_2-stopact_Meal_1)),
         meal2dist = ifelse(sign(meal2dist)== -1,meal2dist+hms::hms(seconds = 0, minutes = 0, hours=24),meal2dist),
         meal2dist = hms::hms(seconds = meal2dist))%>%
  mutate(m2d = as.numeric(meal2dist)/60)%>%
  dplyr::filter(!is.na(m2d))

# specification search - Full Time
meal2_ft1 <- fitdistrplus::fitdist(data = m2$m2d,distr = "norm")
meal2_ft2 <- fitdistrplus::fitdist(m2$m2d,distr = "TF_safe",
                                   start = list(mu=mean(m2$m2d),sigma=sd(m2$m2d),nu=2),
                                   lower=c(100,1,0), upper = c(400,Inf,Inf),method = "mle")
meal2_ft3 <- fitdistrplus::fitdist(m2$m2d,distr = "gompertz",
                                   start = list(shape = 0.01,rate=0.01),
                                   lower=c(0.001,0.001))
meal2_ft4 <- fitdistrplus::fitdist(m2$m2d,distr = "weibull")
meal2_ft5 <- fitdistrplus::fitdist(m2$m2d,distr="gamma")

cdfcomp(list(meal2_ft1,meal2_ft2,meal2_ft3,meal2_ft4,meal2_ft5))
# Gamma fits well.

# specification search - Part Time
meal2_pt1 <- fitdistrplus::fitdist(data=m2pt$m2d,distr = "gamma")
meal2_pt2 <- fitdistrplus::fitdist(m2pt$m2d,distr = "weibull")
# Gamma also fits well.

# Extract Coefficient 
M2params_ptc <- coef(meal2_pt1)
M2params_ftc <- coef(meal2_ft5)

# Sim Check

m2_ptc <- replicate(100,rgamma(500,rate = M2params_ptc[[2]], shape = M2params_ptc[[1]]))
dm2_ptc <- apply(X = m2_ptc,MARGIN = 2,FUN = density,window="e")

meal2denspt <- density(m2pt$m2d)

par(mfrow = c(1, 2),oma = c(0, 1, 3, 0))
plot(meal2denspt,ylim = c(0,.005), 
     main = "Part-Time Students\n(Gamma)",
     bty='n',ylab = "",las=1,cex=.8)
for(i in 1:100){
  lines(dm2_ptc[[i]],col='grey')
}
lines(meal2denspt)

meal2densft <- density(m2$m2d)

m2_ftc <- replicate(100,rgamma(500,rate = M2params_ftc[[2]],shape = M2params_ftc[[1]]))
dm2_ftc <- apply(X = m2_ftc,MARGIN = 2,FUN = density, window="e")

plot(meal2densft,ylim = c(0,.005), 
     main = "Full-Time Students\n(Gamma)",
     bty='n',las=1, ylab = "")
for(i in 1:100){
  lines(dm2_ftc[[i]],col='grey')
}
lines(meal2densft)
mtext("Simulated Density Check: Period Between Meals 1 and 2", side = 3, outer = TRUE, cex = 1.5, font = 2)
mtext("Density", side = 2, outer = TRUE, cex = 1, line = -.1, font = 1)
dev.off()

# ---- Meal 3 ----

m3 <- m_ft %>% 
  mutate(activity = case_when(mnum >= 1 ~ paste("Meal",mnum,sep="_"),
                              snum >= 1 ~ paste("Sleeping",snum,sep="_"), 
                              TRUE ~ activity)) %>%
  dplyr::filter(activity %in% c("Meal_2","Meal_3")) %>%
  dplyr::select(TUCASEID,activity,startact,stopact)%>%
  pivot_wider(id_cols = TUCASEID,names_from = activity,values_from = c(startact,stopact))%>%
  mutate(meal3dist = hms::as_hms(as.numeric(startact_Meal_3-stopact_Meal_2)),
         meal3dist = ifelse(sign(meal3dist)== -1,meal3dist+hms::hms(seconds = 0, minutes = 0, hours=24),
                            meal3dist),
         meal3dist = hms::hms(seconds = meal3dist))%>%
  mutate(m3d = as.numeric(meal3dist)/60)%>%
  dplyr::filter(!is.na(m3d))

m3pt <- m_pt %>% 
  mutate(activity = case_when(mnum >= 1 ~ paste("Meal",mnum,sep="_"),
                              snum >= 1 ~ paste("Sleeping",snum,sep="_"), 
                              TRUE ~ activity)) %>%
  dplyr::filter(activity %in% c("Meal_2","Meal_3")) %>%
  dplyr::select(TUCASEID,activity,startact,stopact)%>%
  pivot_wider(id_cols = TUCASEID,names_from = activity,values_from = c(startact,stopact))%>%
  mutate(meal3dist = hms::as_hms(as.numeric(startact_Meal_3-stopact_Meal_2)),
         meal3dist = ifelse(sign(meal3dist)== -1,meal3dist+hms::hms(seconds = 0, minutes = 0, hours=24),
                            meal3dist),
         meal3dist = hms::hms(seconds = meal3dist))%>%
  mutate(m3d = as.numeric(meal3dist)/60)%>%
  dplyr::filter(!is.na(m3d))

## Very few part time students... use full-time student distribution with censorship rate

# Specification Search
meal3_ft1 <- fitdistrplus::fitdist(data = m3$m3d,distr = "norm")
meal3_ft2 <- fitdistrplus::fitdist(m3$m3d,distr = "TF_safe",
                                   start = list(mu=mean(m3$m3d),sigma=sd(m3$m3d),nu=2),
                                   lower=c(100,1,0), upper = c(400,Inf,Inf),method = "mle")
meal3_ft3 <- fitdistrplus::fitdist(m3$m3d,distr = "weibull")
meal3_ft4 <- fitdistrplus::fitdist(m3$m3d,distr="gamma")

cdfcomp(list(meal3_ft1,meal3_ft2,meal3_ft3,meal3_ft4))
# Normal distribution fits well.

M3params_ftc <- coef(meal3_ft1)

# Sim Check

m3_ftc <- replicate(100,rnorm(500,mean = M3params_ftc[[1]], sd = M3params_ftc[[2]]))
dm3_ftc <- apply(X = m3_ftc,MARGIN = 2,FUN = density,window="e")

meal3densft <- density(m3$m3d)

par(mfrow = c(1, 2),oma = c(0, 1, 3, 0))
plot(meal3densft,ylim = c(0,.005), 
     main = "Full-Time Students\n(Gaussian)",
     bty='n',ylab = "",las=1,cex=.8)
for(i in 1:100){
  lines(dm3_ftc[[i]],col='grey')
}
lines(meal3densft)
mtext("Simulated Density Check: Period Between Meals 2 and 3", side = 3, outer = TRUE, cex = 1.5, font = 2)
mtext("Density", side = 2, outer = TRUE, cex = 1, line = -.1, font = 1)

hist(m3pt$m3d,main = "Part-Time Students \n(Histogram)",xlab= "Minutes")
dev.off()

# ---- All Parameters ----

params <- data.frame(Var = rep(c("Sleep","Meal 1","Meal 2","Meal 3"),each=4),
                     Pop = c(rep(c("Full-Time","Part-Time"), each = 2,lengthout=16)),
                     Params = c(params_ftc,params_ptc,Mparams_ftc,Mparams_ptc,M2params_ftc,M2params_ptc,M3params_ftc,NA,NA),
                     Labs = c("mean","sd","shape","scale",rep(c("rate","pi"),2),rep(c("shape","rate"),2),"mean","sd","NA","NA"),
                     Dist = c("Gaussian","Gaussian","Weibull","Weibull",
                              rep("Zero-Inflated Exponential",4),
                              rep("Gamma",4),
                              rep(c("Gaussian","NA"),each=2)))%>%
  arrange(Pop,Var)%>%
  dplyr::select(Pop, Var, Dist, Labs, Params)%>%
  rbind(data.frame(Pop = c("Full-Time","Part-Time"),
        Var = c("Meal 3", "Meal 3"),
        Dist = NA,
        Labs = rep("Censorship Rate",2),
        Params = c(round((207-68)/207,3),round((53-21)/53,2))))%>%
  arrange(Pop,Var)

write_csv(params,here("data/proc_data","02_SleepMealParams.csv"))

# ---- Activity Distributions ----

# Total minutes recorded per person
pptotals <- diary_ftc %>% 
  summarize(totmin = sum(TUACTDUR24),.by = TUCASEID)
pptotals_pt <- diary_ptc %>% 
  summarize(totmin = sum(TUACTDUR24),.by = TUCASEID)

# Observed distribution, complete cases with 0's for unobserved activities
# Full time students
actdist <- diary_ftc %>%
  summarize(mmin = sum(TUACTDUR24),
            .by = c(TUCASEID,activity))%>%
  complete(activity, TUCASEID, fill = list(mmin = 0))%>%
  left_join(pptotals)%>%
  mutate(props = mmin/totmin)%>%
  summarize(mprop_ft = mean(props),.by=activity)

# Part time students 
actdistpt <- diary_ptc %>%
  summarize(mmin = sum(TUACTDUR24),
            .by = c(TUCASEID,activity))%>%
  complete(activity, TUCASEID, fill = list(mmin = 0))%>%
  left_join(pptotals_pt)%>%
  summarize(mprop_pt = sum(mmin)/sum(totmin),.by=activity)

# Combine 
actd <- actdist %>%
  full_join(actdistpt)%>%
  mutate(mprop_pt = ifelse(is.na(mprop_pt),0,mprop_pt),
         mprop_ft = round(mprop_ft,4),
         mprop_pt = round(mprop_pt,4))

write_csv(actd,here("data/proc_data","02_ActDist.csv"))


actd <- actd %>%
  mutate(diff = mprop_ft-mprop_pt,
         odds = mprop_ft/mprop_pt)

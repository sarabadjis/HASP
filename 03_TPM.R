## Transition Matrix
# August 2026
# Sophia Arabadjis

library(here)
library(tidyverse)
# Call instead of load
#library(fitdistrplus)
#library(hms)

source(here("code","00_Source.R"))

diary_ftc <- readRDS(here("data/proc_data","01_FTS.rds"))
diary_ptc <- readRDS(here("data/proc_data","01_PTS.rds"))
sdist <- read_csv(here("data/proc_data","02_StatDist.csv"))

actlist <- unique(unique(diary_ftc$activity),unique(diary_ptc$activity))%>% sort()

# Intialize
ftc_id <- unique(diary_ftc$TUCASEID)
ptc_id <- unique(diary_ptc$TUCASEID)

# TMPs
tmp_ftc <- tm_tab.f(df = diary_ftc,id_list=ftc_id,acts=actlist,verbose=T)
tmp_ptc <- tm_tab.f(df = diary_ptc,id_list=ptc_id,acts=actlist,verbose=T)

plot(tmp_ftc$Trace[2:length(tmp_ftc$Trace)],type = 'l',col='skyblue', 
     ylab = "Tolerance (average difference)",
     xlab = "Iterations",
     main = "Transition Probability Matrix Trace",
     ylim = c(0,0.003))
lines(tmp_ptc$Trace[2:length(tmp_ptc$Trace)], col = 'orange')
abline(h = 0.0001,lty=2,col='black')

## ---- Solve for Stationary Distributions ----
ftc_eigen <- eigen(t(tmp_ftc$Product))
ftc_vec <- Re(ftc_eigen$vectors[, 1])
ftc_statdist <- ftc_vec/sum(ftc_vec)

ptc_eigen <- eigen(t(tmp_ptc$Product))
ptc_vec <- Re(ptc_eigen$vectors[, 1])
ptc_statdist <- ptc_vec/sum(ptc_vec)

statdist <- t(rbind(ftc_statdist%*%tmp_ftc$Raw,
                  ptc_statdist%*%tmp_ptc$Raw))%>%
  as.data.frame()

write_csv(statdist,here("Output","03_StationaryDistributions.csv"))

## ---- Multinomial Sampling ----

sample_MN <- MNsim_sampler_P.f(nsize = 30,
                               tsize=3, 
                               tvec=c(2*6,24,2*18,48),
                               u = 48,
                               probvector = statdist$V1,
                               activityvector = rownames(statdist),
                               output = "samplematrix")
colrs <- c("goldenrod","goldenrod4",
           "skyblue","skyblue4",
           "thistle","thistle4",
           "wheat","wheat4",
           "tomato4","tomato",
           "seagreen1","seagreen4",
           'olivedrab','olivedrab3')

named_colors <- setNames(as.character(colrs), as.character(actlist))


data.frame(sample_MN) %>%
  pivot_longer(cols = starts_with("X"),
               names_to = "Simulated Person",
               values_to = "Activity") %>%
  #mutate(Activity = factor(Activity,labels=actlist,levels=1:14))%>%
  arrange(`Simulated Person`) %>%
  mutate(Time = rep(48:1,30),
         Person = rep(1:30,each=48)) %>%
  ggplot(aes(Time,Person,fill = Activity))+
  geom_tile()+
  scale_fill_manual(values=named_colors,drop=F)+
  theme_classic()+
  scale_y_continuous(name = "Simulated Person", breaks=0:30)+
  scale_x_continuous(name = "Time",breaks = c(0,12,24,36,48))+
  theme(axis.text = element_blank(),
        plot.title = element_text(hjust=.5))+
  ggtitle("A: Multinomial Generated Activities\n(Full-Time Students)")

## ---- Monte Carlo Sampling ----

# Rearrange to have sleep first row and first column
mmsleep <- tmp_ftc$Product[c(9,1:8,10:14),c(9,1:8,10:14)]
avsleep <- paste(actlist[c(9,1:8,10:14)])

set.seed(8282025)

sample_MC <- MN_corr_sampler2.f(mm = mmsleep, 
                               av = avsleep,
                               nsize=30,
                               tsize = 3,
                               tvec = c(15,30,48),
                               reps = 48,
                               output = "samplematrix")

data.frame(sample_MC) %>%
  pivot_longer(cols = starts_with("X"),
               names_to = "Simulated Person",
               values_to = "Activity") %>%
  arrange(`Simulated Person`) %>%
  mutate(Time = rep(49:1,30),
         Person = rep(1:30,each=49)) %>%
  ggplot(aes(Time,Person,fill = Activity))+
  geom_tile(show.legend = F)+
  scale_fill_manual(values=named_colors)+
  theme_classic()+
  scale_y_continuous(name = "Simulated Person", breaks=0:30)+
  scale_x_continuous(name = "Time",breaks = c(0,12,24,36,48))+
  theme(axis.text = element_blank(),
        plot.title = element_text(hjust=.5))+
  ggtitle("B: Monte Carlo Generated Activities\n(Full-Time Students)")

## ---- Period-Based Transition Probabilities ----

diary_ftc2 <- diary_ftc %>%
  group_by(TUCASEID,activity)%>%
  arrange(TUACTIVITY_N)%>%
  mutate(countnum = row_number())%>%
  ungroup()%>%
  mutate(act = activity,
         activity = as.factor(ifelse(act=="Eating and Drinking",
                           paste("Eating and Drinking",countnum,sep=" "),paste(act))))

diary_ptc2 <- diary_ptc %>%
  group_by(TUCASEID,activity)%>%
  arrange(TUACTIVITY_N)%>%
  mutate(countnum = row_number())%>%
  ungroup()%>%
  mutate(act = activity,
         activity = as.factor(ifelse(act=="Eating and Drinking",
                      paste("Eating and Drinking",countnum,sep=" "),paste(act))))

al <- unique(unique(diary_ftc2$activity),unique(diary_ptc2$activity))%>% sort()

write_csv(data.frame(activity = al),here("data/proc_data","03_ActivityList.csv"))

# Period Diaries
ft_period <- period_diaries.f(diary_ftc2)
pt_period <- period_diaries.f(diary_ptc2)

alsmall <- paste(al[-c(1:3,11)])

# Tabulate

FTP1 <- tm_Period_tab.f(df = ft_period$Period1,
                        id_list=unique(ft_period$Period1$TUCASEID),
                        acts = alsmall,
                        Period = 1)
FTP2 <- tm_Period_tab.f(df = ft_period$Period2,
                        id_list=unique(ft_period$Period2$TUCASEID),
                        acts = alsmall,
                        Period = 2)
FTP3 <- tm_Period_tab.f(df = ft_period$Period3,
                        id_list=unique(ft_period$Period3$TUCASEID),
                        acts = alsmall,
                        Period = 3)
PTP1 <- tm_Period_tab.f(df = pt_period$Period1,
                        id_list=unique(pt_period$Period1$TUCASEID),
                        acts = alsmall,
                        Period = 1)
PTP2 <- tm_Period_tab.f(df = pt_period$Period2,
                        id_list=unique(pt_period$Period2$TUCASEID),
                        acts = alsmall,
                        Period = 2)
PTP3 <- tm_Period_tab.f(df = pt_period$Period3,
                        id_list=unique(pt_period$Period3$TUCASEID),
                        acts = alsmall,
                        Period = 3)


write_rds(list(FTP1 = FTP1, PTP1 = PTP1,
               FTP2 = FTP2, PTP2 = PTP2,
               FTP3 = FTP3, PTP3 = PTP3),
          here("Output","03_TPM.rds"))


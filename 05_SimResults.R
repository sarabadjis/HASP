## Simulations
# September 2026
# Sophia Arabadjis

library(here)
library(tidyverse)
library(simsurv)
library(flexsurv)
library(patchwork)

source(here("code","00_Source.R"))

FT_sim <- readRDS(here("Output","04_FuTiSim.rds"))
PT_sim <- readRDS(here("Output","04_PaTiSim.rds"))

# ----------------
# Data Cleaning
# ----------------

# Cut data across day lines (multiples of 1440):
FTk <- daysect.f(FT_sim)
PTk <- daysect.f(PT_sim)

FTk <- FTk %>% arrange(id,tstart)
PTk <- PTk %>% arrange(id,tstart)

## Get minute-by-minute activity breakdowns:

g1 <- FTk %>%
  dplyr::select(id,tstart,tstop,activity)%>%
  mutate(tstart = round(tstart,0),
         tstop = round(tstop,0),
         day = floor((tstart+1)/1440))

g2 <- PTk %>%
  dplyr::select(id,tstart,tstop,activity)%>%
  mutate(tstart = round(tstart,0),
         tstop = round(tstop,0),
         day = floor((tstart+1)/1440))

## Expand Grid to get full minute-by-minute break down for each group
# Full-Time 1:
daychoice <- expand.grid(tstart = seq(2880,2880+1440,1),
                         id = seq(1,1000,1),
                         day = floor(2880/1440),
                         activity = NA)

daychoice <- daychoice %>%
  left_join(
    g1 %>% select(id, tstart_g = tstart, tstop_g = tstop, day_g = day, activity),
    by = join_by(id, tstart >= tstart_g, tstart <= tstop_g, day == day_g)
  ) %>%
  mutate(activity = coalesce(activity.y, activity.x)) %>%
  select(-activity.x, -activity.y, -tstart_g, -tstop_g)%>%
  group_by(id,tstart)%>% # Force transtion at second minute
  mutate(g = row_number())%>%
  ungroup()%>%
  dplyr::filter(g==1)

# Part Time 
daychoice2 <- expand.grid(tstart = seq(2880,2880+1440,1),
                          id = seq(1,1000,1),
                          day = floor(2880/1440),
                          activity = NA)

daychoice2 <- daychoice2 %>%
  left_join(
    g2 %>% select(id, tstart_g = tstart, tstop_g = tstop, day_g = day, activity),
    by = join_by(id, tstart >= tstart_g, tstart <= tstop_g, day == day_g)
  ) %>%
  mutate(activity = coalesce(activity.y, activity.x)) %>%
  select(-activity.x, -activity.y, -tstart_g, -tstop_g)%>%
  group_by(id,tstart)%>%
  mutate(g = row_number())%>%
  ungroup()%>%
  dplyr::filter(g==1)

saveRDS(daychoice2,here("Output","daychoicePT.rdata"))
saveRDS(daychoice,here("Output","daychoiceFT.rdata"))

## Summarize:
mls <- c(paste("Eating and Drinking",1:3,sep=" "))
# Group 1
daychoice_sum <- daychoice %>%
  mutate(activity = ifelse(activity %in% mls,"Eating and Drinking",activity))%>%
  summarize(n = n()/1000, .by = c(tstart,activity)) %>% 
  arrange(tstart,activity)

# Group 2
daychoice_sum2 <- daychoice2 %>%
  mutate(activity = ifelse(activity %in% mls,"Eating and Drinking",activity))%>%
  summarize(n =n()/1000,.by = c(tstart,activity)) %>% 
  arrange(tstart,activity)

# --------
# Visuals 
# -------

colrs <- c("goldenrod","goldenrod4",
           "skyblue","skyblue4",
           "thistle","thistle4",
           "wheat","wheat4",
           "tomato4","tomato",
           "seagreen4","seagreen1",
           'olivedrab','olivedrab3')

named_colors <- setNames(as.character(colrs), unique(daychoice_sum$activity)[1:14] %>% sort())

# Generate Activity Profile Figure for Group 1

fig1 <- daychoice_sum %>%
  group_by(tstart)%>%
  arrange(tstart,activity)%>%
  mutate(yend = cumsum(n),
         ystart = yend-diff(c(0,yend),lag = 1)) %>%
  dplyr::filter(!is.na(activity))%>% #First row
  ggplot()+
  geom_segment(aes(x=tstart,xend=tstart,
                   y = ystart, yend=yend, 
                   colour = activity))+
  scale_color_manual(values=named_colors,name="Activity")+
  scale_x_continuous(breaks = c(2880,3240,3600,3960,4320), labels = c(0,6,12,18,24))+
  labs(x = "Time (hours)", y = expression(P[j](t)), title = "Full-Time Student")+
  theme_classic()

# Generate Activity Profile Figure for Group 2

fig2 <-  daychoice_sum2 %>%
  group_by(tstart)%>%
  arrange(tstart,activity)%>%
  mutate(yend = cumsum(n),
         ystart = yend-diff(c(0,yend),lag = 1)) %>%
  dplyr::filter(!is.na(activity))%>%
  ggplot()+
  geom_segment(aes(x=tstart,xend=tstart,
                   y = ystart, yend=yend, 
                   colour = activity),show.legend = F)+
  scale_color_manual(values=named_colors)+
  scale_x_continuous(breaks = c(2880,3240,3600,3960,4320), labels = c(0,6,12,18,24))+
  labs(x = "Time (hours)", y = expression(P[j](t)), title = "Part-Time Student")+
  theme_classic()


## ---- Survival Curves ----
## Check Sleep

group1_sleep <- FT_sim %>%
  dplyr::filter(activity=="Sleeping")%>%
  mutate(diff_t = tstop-tstart,
         diff_s = 0,
         event = 1,
         group="Full-Time")%>%
  dplyr::select(id,diff_s,diff_t,event,group)
group2_sleep <- PT_sim %>%
  dplyr::filter(activity=="Sleeping")%>%
  mutate(diff_t = tstop-tstart,
         diff_s = 0,
         event = 1,
         group="Part-Time")%>%
  dplyr::select(id,diff_s,diff_t,event,group)

sleep <- rbind(group1_sleep, group2_sleep)

fitsleepgroup <-survfit(Surv(sleep$diff_s/60, sleep$diff_t/60,event=sleep$event)~sleep$group)
fitsleeptotal <- survfit(Surv(sleep$diff_s/60, sleep$diff_t/60,event=sleep$event)~1)

plot(fitsleeptotal, 
     main = "Sleep Durations",xlab = "Hours", font.main=1, xlim=c(0,15),
     col = c("tomato4",alpha("tomato4",.65),alpha("tomato4",.65)))
for (i in 1:100){
  i <- i
  s <- sample(1:dim(sleep)[1],100,replace=F)
  tmp1 <- group1_sleep[s,]
  tmp2 <- group2_sleep[s,]
  fit1 <- survfit(Surv(tmp1$diff_s/60, tmp1$diff_t/60,
                      event=tmp1$event)~ 1)
  fit2 <- survfit(Surv(tmp2$diff_s/60, tmp2$diff_t/60,
                       event=tmp2$event)~ 1)
  lines(fit1[1],col = "skyblue1",add=T,conf.int = F, lwd=.4)
  lines(fit2[1],col = "orange1",add=T,conf.int = F, lwd=.4)
}
lines(fitsleeptotal, 
      col = "tomato4", add=T)
lines(fitsleepgroup,
      col = c("skyblue4","orange3"), add=T)
segments(x0 = 0, x1 = 9.04,y0=0.5,y1=0.5, col = "black", lwd=1,lty=2)
segments(x0 = 9.34, x1 = 9.07,y0=0,y1=0.5, col = "black", lwd=1,lty=2)
text(x=13,y=.85,paste("Median Sleep Time:",round(median(fitsleeptotal)[1],1),"hours"))
text(x=13,y=.8,paste("Full-Time Students:",round(median(fitsleepgroup)[1],1),"hours"))
text(x=13,y=.75,paste("Part-Time Students:",round(median(fitsleepgroup[2]),1),"hours"))
legend("bottomleft", 
       legend=c("All Students","Full-Time Students","Part-Time Students",
                "Random Full-Time Student Samples (N=100)",
                "Random Part-Time Student Samples (N=100)"),
       col=c("tomato4",'skyblue4','orange3','skyblue1','orange1'),lwd=c(1,1,1,.5,.5), 
       cex=.75, bty='n',lty=1)

## Meals

group1_meals <- FT_sim %>%
  dplyr::filter(str_detect(activity,"Eating"))%>%
  mutate(diff_t = tstart %% 1440 / 60,
         diff_s = 0,
         event = 1)%>%
  dplyr::select(id,diff_s,diff_t,event,activity)
group2_meals <- PT_sim %>%
  dplyr::filter(str_detect(activity,"Eating"))%>%
  mutate(diff_t = tstart %% 1440 / 60,
         diff_s = 0,
         event = 1)%>%
  dplyr::select(id,diff_s,diff_t,event,activity)


a <- survfit(Surv(group1_meals$diff_s,group1_meals$diff_t, event=group1_meals$event)~group1_meals$activity)
b <- survfit(Surv(group2_meals$diff_s,group2_meals$diff_t, event=group2_meals$event)~group2_meals$activity)

plot(a,
     col = c("#D95F02","orange", "tomato3"),conf.int = 0.95, xlab = "Hours",
     main = "Meal Occurance Full-Time Students", font.main= 1)
for (i in 1:50){
  i <- i
  s <- sample(1:dim(group1_meals)[1],100,replace=F)
  tmp <- group1_meals[s,]
  fit <- survfit(Surv(tmp$diff_s, tmp$diff_t,
                      event=tmp$event)~ tmp$activity)
  lines(fit[1],col = "lightgrey",add=T,conf.int = F, lwd=.4)
  lines(fit[2],col = "lightgrey",add=T,conf.int = F, lwd=.4)
  lines(fit[3],col = "lightgrey",add=T,conf.int = F,lwd=.4)
}
lines(a[1:3],col = c("#D95F02","orange", "tomato3"))
segments(x0 = 0, x1 = median(a)[3],y0=0.5,y1=0.5, col = "tomato3", lwd=2)
segments(x0 = 0, x1 = median(a)[2],y0=0.5,y1=0.5, col = "orange", lwd=2)
segments(x0 = 0, x1 = median(a)[1],y0=0.5,y1=0.5, col = "#D95F02", lwd=2)
text(x = 18.5, y =0.55,round(median(a)[3],1),col="tomato3")
text(x = 13,y=0.55,round(median(a)[2],1),col="orange")
text(x = 8,y= 0.55,round(median(a)[1],1),col="#D95F02")
legend("bottomleft",legend=c("First","Second","Third","Random Samples (N=100)"), lty = 1, 
       col = c("#D95F02","orange","tomato3","lightgrey"),bty='n',
       cex=.7, lwd=c(2,2,2,1))

plot(b,
     col = c("#D95F02","orange", "tomato3"),conf.int = 0.95, xlab = "Hours",
     main = "Meal Occurance Part-Time Students", font.main= 1)
for (i in 1:50){
  i <- i
  s <- sample(1:dim(group2_meals)[1],100,replace=F)
  tmp <- group2_meals[s,]
  fit <- survfit(Surv(tmp$diff_s, tmp$diff_t,
                      event=tmp$event)~ tmp$activity)
  lines(fit[1],col = "lightgrey",add=T,conf.int = F, lwd=.4)
  lines(fit[2],col = "lightgrey",add=T,conf.int = F, lwd=.4)
  lines(fit[3],col = "lightgrey",add=T,conf.int = F,lwd=.4)
}
lines(b[1:3],col = c("#D95F02","orange", "tomato3"))
segments(x0 = 0, x1 = median(b)[3],y0=0.5,y1=0.5, col = "tomato3", lwd=2)
segments(x0 = 0, x1 = median(b)[2],y0=0.5,y1=0.5, col = "orange", lwd=2)
segments(x0 = 0, x1 = median(b)[1],y0=0.5,y1=0.5, col = "#D95F02", lwd=2)
text(x = 18, y =0.55,round(median(b)[3],1),col="tomato3")
text(x = 13,y=0.55,round(median(b)[2],1),col="orange")
text(x = 7.5,y= 0.55,round(median(b)[1],1),col="#D95F02")
legend("bottomleft",legend=c("First","Second","Third","Random Samples (N=100)"), lty = 1, 
       col = c("#D95F02","orange","tomato3","lightgrey"),bty='n',
       cex=.7, lwd=c(2,2,2,1))

## Household Differences
group1_sleep_start <- FT_sim %>%
  dplyr::filter(activity=="Sleeping")%>%
  dplyr::filter(tstart > 1440)%>%
  dplyr::filter(tstart < 2880+360)%>%
  mutate(diff_t = tstart %% 1440,
         event = 1,
         group = 1)%>%
  group_by(id)%>%
  slice_max(tstart)%>%
  dplyr::select(id,diff_t,HH, group)
group2_sleep_start <- PT_sim %>%
  dplyr::filter(activity=="Sleeping")%>%
  dplyr::filter(tstart > 1440)%>%
  dplyr::filter(tstart < 2880+360)%>%
  mutate(diff_t = tstart %% 1440,
         diff_s = 0,
         event = 1,
         group = 2,
         HH = HH+200)%>%
  group_by(id)%>%
  slice_max(tstart)%>%
  dplyr::select(id,diff_t,HH, group)

sleep_var <- rbind(group1_sleep_start, group2_sleep_start)%>%
  mutate(HH = factor(HH),
         diff_c = ifelse(diff_t<400,diff_t+1440,diff_t))

c <- lm(diff_c ~ as.factor(group), data = sleep_var)
d <- lm(diff_c ~ as.factor(group) + as.factor(HH), data = sleep_var)
e <- lme4::lmer(diff_c ~ 1 + (1|HH), data= sleep_var %>% mutate(HH = as.factor(HH),
                                                                  group = as.factor(group)))
anova_mods <- anova(c,d)

resdat <- data.frame(Residuals = c(c$residuals,d$residuals,residuals(e)),
                     Fitted = c(c$fitted.values,d$fitted.values,fitted.values(e)),
                     Type = rep(c("A) Reduced Model","B) Fixed Effect Model","C) Random Effect Model"),each=2000),
                     Index = rep(1:2000,3))

colchoice <- c("A) Reduced Model" = "forestgreen",
               "B) Fixed Effect Model" = 'wheat4',
               "C) Random Effect Model" = 'plum4')

figres <- resdat %>%
  ggplot(aes(Index,Residuals,group = Type))+
  geom_point(aes(shape = Type, col = Type),size=.5)+
  scale_color_manual(values= colchoice)+
  geom_smooth(aes(Index,Residuals),method = loess,col = "black")+
  facet_wrap(~Type)+
  theme_linedraw()+
  theme(legend.position = 'none')

## Table of Observed Differences:

FT_diff <- g1 %>%
  dplyr::filter(day==2)%>%
  mutate(dur = tstop-tstart,
         a = ifelse(str_detect(activity,"Eating"),"Eating and Drinking",activity))%>%
  summarize(ndur = sum(dur),.by=c(id,a))%>%
  mutate(nprop = ndur/(1440))%>%
  summarize(npropm = sum(nprop)/1000,.by=a)
PT_diff <- g2 %>%
  dplyr::filter(day==2)%>%
  mutate(dur = tstop-tstart,
         a = ifelse(str_detect(activity,"Eating"),"Eating and Drinking",activity))%>%
  summarize(ndur = sum(dur),.by=c(id,a))%>%
  mutate(nprop = ndur/(1440))%>%
  summarize(npropm = sum(nprop)/1000,.by=a)

Tab_diff <- FT_diff %>%
  full_join(PT_diff,by = "a",suffix = c(".ft",".pt"))%>%
  mutate(npropm.pt = ifelse(is.na(npropm.pt),0,npropm.pt),
         ObDiff = abs(npropm.pt-npropm.ft),
         OR = npropm.ft/npropm.pt)%>%
  arrange(a)

sum(round(Tab_diff$npropm.ft,3))
sum(round(Tab_diff$npropm.pt,3))

## ---- Check of Individual Activity Realism:

# Sample 10 from each group

ftss <- sample(1:1000,10)
ptss <- sample(1:1000,10)

FTSS <- g1[g1$id %in% ftss,]
PTSS <- g2[g2$id %in% ptss,]

# Plot
tots <- FTSS %>%
  dplyr::filter(day==2)%>%
  group_by(id)%>%
  mutate(rn = cur_group_id())%>%
  ungroup()%>%
  mutate(activity = ifelse(str_detect(activity,"Eating"),
                           "Eating and Drinking",
                           activity),
         group="Full-Time Students")%>%
  arrange(id,tstart)%>%
  rbind(PTSS %>%
          dplyr::filter(day==2)%>%
          group_by(id)%>%
          mutate(rn = cur_group_id())%>%
          ungroup()%>%
          mutate(activity = ifelse(str_detect(activity,"Eating"),
                                   "Eating and Drinking",
                                   activity),
                 group = "Part-Time Students")%>%
          arrange(id,tstart))%>%
  mutate(activity = factor(activity,levels=unique(daychoice_sum$activity)[1:14] %>% sort(),
                    labels = unique(daychoice_sum$activity)[1:14] %>% sort()))


fcheck <- tots %>%
  dplyr::filter(group=="Full-Time Students")%>%
  ggplot()+
  geom_segment(aes(x=tstart,xend=tstop,y=rn,yend=rn,col=activity),lwd=4,show.legend = T) +
  scale_color_manual(values = named_colors, name="Activity", drop=FALSE)+
  scale_y_continuous(breaks=c(1:10),name="Student No.")+
  scale_x_continuous(breaks = seq(2880,4320,length.out=7),
                     labels = c(0,4,8,12,16,20,24),
                     name = "Time")+
  ggtitle("Full-Time Students")+
  theme_classic()

pcheck <- tots %>%
  dplyr::filter(group=="Part-Time Students")%>%
  ggplot()+
  geom_segment(aes(x=tstart,xend=tstop,y=rn,yend=rn,col=activity),lwd=4,show.legend = F) +
  scale_color_manual(values = named_colors,name="Activity")+
  scale_y_continuous(breaks=c(1:10),name="Student No.")+
  scale_x_continuous(breaks = seq(2880,4320,length.out=7),
                     labels = c(0,4,8,12,16,20,24),
                     name = "Time")+
  ggtitle("Part-Time Students")+
  theme_classic()

fcheck + pcheck + plot_annotation(title = "10 Simulated Individual Activity Schedules")

## Person 6 check

F6 <- tots[tots$rn == 6 & tots$group == "Full-Time Students",]
P6 <- tots[tots$rn == 6 & tots$group == "Part-Time Students",]

F6 <- F6 %>%
  dplyr::filter(day==2)%>%
  group_by(id, activity, group_id = consecutive_id(activity))%>%
  summarize(tstop = max(tstop),
            tstart = min(tstart))%>%
  arrange(tstart)%>%
  mutate(tstart = tstart %% 1440,
         tstop = tstop %% 1440,
         dur = tstop-tstart)%>%
  ungroup()%>%
  dplyr::select(activity,tstart,tstop,dur)

P6 <- P6 %>%
  dplyr::filter(day==2)%>%
  group_by(id, activity, group_id = consecutive_id(activity))%>%
  summarize(tstop = max(tstop),
            tstart = min(tstart))%>%
  arrange(tstart)%>%
  mutate(tstart = tstart %% 1440,
         tstop = tstop %% 1440,
         dur = tstop-tstart)%>%
  ungroup()%>%
  dplyr::select(activity,tstart,tstop,dur)

ttt <- cbind(F6[1:10,],P6[1:10,])

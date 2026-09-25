## Comparison to Observed Data
# September 2026
# Sophia Arabadjis

library(here)
library(tidyverse)
library(simsurv)
library(flexsurv)
library(patchwork)

source(here("code","00_Source.R"))


## Read in Data 

fts <- readRDS(here("data/proc_data","01_FTS.rds"))
pts <- readRDS(here("data/proc_data","01_PTS.rds"))

colrs <- c("goldenrod","goldenrod4",
           "skyblue","skyblue4",
           "thistle","thistle4",
           "wheat","wheat4",
           "tomato4","tomato",
           "seagreen1","seagreen4",
           'olivedrab','olivedrab3')

named_colors <- setNames(as.character(colrs), 
                         unique(daychoice_sum$activity)[1:14] %>% sort())

## Clean Data for easier expansion: 

fts <- fts %>%
  mutate(startactm = hms::hms(seconds = as.numeric(startact)),
         stopactm = startactm+hms::hms(minutes = as.numeric(TUACTDUR24)),
         stopactm = hms::hms(seconds=as.numeric(stopactm)))%>%
  rename(id = TUCASEID)

pts <- pts %>%
  mutate(startactm = hms::hms(seconds = as.numeric(startact)),
         stopactm = startactm+hms::hms(minutes = as.numeric(TUACTDUR24)),
         stopactm = hms::hms(seconds=as.numeric(stopactm)))%>%
  rename(id = TUCASEID)

## Create minute-by-minute
daymin <- expand.grid(tstart = seq(240,240+1440,1),
                          id = unique(fts$id),
                          activity = NA)%>%
  mutate(tstart = hms::hms(minutes = as.numeric(tstart)))

daymin2 <- expand.grid(tstart = seq(240,240+1440,1),
                      id = unique(pts$id),
                      activity = NA)%>%
  mutate(tstart = hms::hms(minutes = as.numeric(tstart)))

## Joins
dayminfts <- daymin %>%
  left_join(
    fts %>% select(id,startactm, stopactm, activity),
    by = join_by(id, tstart >= startactm, tstart <= stopactm)
  ) %>%
  mutate(activity = coalesce(activity.y, activity.x)) %>%
  select(-activity.x, -activity.y, -startactm, -stopactm)%>%
  group_by(id,tstart)%>%
  mutate(g = row_number())%>%
  ungroup()%>%
  dplyr::filter(g==1)

dayminpts <- daymin2 %>%
  left_join(
    pts %>% select(id,startactm, stopactm, activity),
    by = join_by(id, tstart >= startactm, tstart <= stopactm)
  ) %>%
  mutate(activity = coalesce(activity.y, activity.x)) %>%
  select(-activity.x, -activity.y, -startactm, -stopactm)%>%
  group_by(id,tstart)%>%
  mutate(g = row_number())%>%
  ungroup()%>%
  dplyr::filter(g==1)

## Normalize 

dayminfts_tots <- dayminfts %>%
  dplyr::filter(!is.na(activity))%>%
  group_by(tstart)%>%
  summarize(n_tot = n())

dayminpts_tots <- dayminpts %>%
  dplyr::filter(!is.na(activity))%>%
  group_by(tstart)%>%
  summarize(n_tot = n())

dayminfts_sum <- dayminfts %>%
  summarize(n = n(), .by = c(tstart,activity)) %>% 
  full_join(dayminfts_tots)%>%
  mutate(n = n/n_tot)%>%
  arrange(tstart,activity)

dayminpts_sum <- dayminpts %>%
  summarize(n = n(), .by = c(tstart,activity)) %>% 
  full_join(dayminpts_tots)%>%
  mutate(n = n/n_tot)%>%
  arrange(tstart,activity)

# Plot

fig1 <- dayminfts_sum %>%
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
  scale_x_continuous(breaks = hms::hms(hours = seq(4,28,4)), labels = hms::hms(hours=seq(4,28,4)))+
  labs(x = "Time (hours)", y = expression(P[j](t)), title = "ATUS Full-Time Students (n=207)")+
  theme_classic()

fig2 <- dayminpts_sum %>%
  group_by(tstart)%>%
  arrange(tstart,activity)%>%
  mutate(yend = cumsum(n),
         ystart = yend-diff(c(0,yend),lag = 1)) %>%
  dplyr::filter(!is.na(activity))%>% #First row
  ggplot()+
  geom_segment(aes(x=tstart,xend=tstart,
                   y = ystart, yend=yend, 
                   colour = activity),show.legend = F)+
  scale_color_manual(values=named_colors,name="Activity")+
  scale_x_continuous(breaks = hms::hms(hours = seq(4,28,4)), labels = hms::hms(hours=seq(4,28,4)))+
  labs(x = "Time (hours)", y = expression(P[j](t)), title = "ATUS Part-Time Students (n=53)")+
  theme_classic()

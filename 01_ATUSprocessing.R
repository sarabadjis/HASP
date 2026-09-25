## Processing ATUS Data
# August 2026
# Sophia Arabadjis

library(here)
library(tidyverse)

# ---- Read in Data ----

df <- read.table(file = here("data/atusact-2024","atusact_2024.dat"), header = T, sep = ",", colClasses = "character")
dfrost <- read.table(file = here("data/atusrost-2024","atusrost_2024.dat"), header = T, sep = ",", colClasses = "character")
dfresp <- read.table(file = here("data/atusresp-2024","atusresp_2024.dat"), header = T, sep = ",", colClasses = "character")

# ---- Clean and Filter to Full time and part time students ----

dfrost <- dfrost %>%
  mutate(TEAGE = as.numeric(TEAGE),
         TESEX = as.factor(TESEX))

dfstud <- dfresp %>%
  dplyr::filter(TESCHENR == 1)%>%
  dplyr::select(-TULINENO)

dfstud_details <- dfrost %>%
  dplyr::filter(TUCASEID %in% dfstud$TUCASEID & (TEAGE <= 23 | TEAGE >=17))%>%
  left_join(dfstud)

# Full time, Part time
ft_stud <- dfstud_details %>%
  dplyr::filter(TESCHFT == 1)
pt_stud <- dfstud_details %>%
  dplyr::filter(TESCHFT == 2)

# Create Diaries of full time and part time students 
diary_ft <- df %>%
  dplyr::filter(TUCASEID %in% ft_stud$TUCASEID)%>%
  dplyr::select(c(TUCASEID,TUACTIVITY_N,TRCODE,TUACTDUR24,TUSTARTTIM,TUSTOPTIME))
diary_pt <- df %>%
  dplyr::filter(TUCASEID %in% pt_stud$TUCASEID)%>%
  dplyr::select(c(TUCASEID,TUACTIVITY_N,TRCODE,TUACTDUR24,TUSTARTTIM,TUSTOPTIME))

#### ---- Recode into Useful Categories -----
# Remove:
# "helping non-household members"
# "services (legal, etc.)",
# "household services (e.g. gardener)",
# "government & civic service",
# "volunteer activities"
# "unable to code"
# Code as "Other"

diary_ftc <- diary_ft %>%
  mutate(activity = case_when(str_detect(TRCODE,"^0101")~ "Sleeping",
                              str_detect(TRCODE,"^0102") ~ "Personal Activities",
                              str_detect(TRCODE,"^0103") ~ "Personal Activities",
                              str_detect(TRCODE,"^0104") ~ "Personal Activities",
                              str_detect(TRCODE,"^0105") ~ "Personal Activities",
                              str_detect(TRCODE,"^0199") ~ "Personal Activities",
                              str_detect(TRCODE,"^02") ~ "Household Activities",
                              str_detect(TRCODE,"^03") ~ "Helping HH Members",
                              str_detect(TRCODE,"^05") ~ "Work",
                              str_detect(TRCODE,"^06") ~ "School",
                              str_detect(TRCODE,"^07") ~ "Shopping",
                              str_detect(TRCODE,"^11") ~ "Eating and Drinking",
                              str_detect(TRCODE,"^12") ~ "Socializing and Leisure",
                              str_detect(TRCODE,"^13") ~ "Sports, Exercise, Recreation",
                              str_detect(TRCODE,"^14") ~ "Religious/Spiritual Activities",
                              str_detect(TRCODE,"^16") ~ "Telephone Calls",
                              str_detect(TRCODE,"^18") ~ "Traveling/Commuting",
                              TRUE ~ "Other"))%>%
  mutate(TUACTDUR24 = as.numeric(TUACTDUR24),
         startact = hms::as_hms(TUSTARTTIM),
         stopact = startact,
         stopact = stopact + hms::hms(minutes = TUACTDUR24),
         stopact = hms::hms(seconds = as.numeric(stopact)),
         stopact = ifelse(stopact > hms::as_hms("24:00:00"),
                          as.numeric(stopact) - (24*60*60),
                          as.numeric(stopact)),
         stopact = hms::hms(seconds = stopact),
         activity = as.factor(activity))%>%
  group_by(TUCASEID, activity, group_id = consecutive_id(activity)) %>%
  summarise(TUACTDUR24 = sum(TUACTDUR24),
            startact = min(startact),
            stopact = max(stopact),
            TUACTIVITY_N = mean(as.numeric(TUACTIVITY_N)), .groups = "drop") %>%
  arrange(group_id) %>%
  select(-group_id)

diary_ptc <- diary_pt %>%
  mutate(activity = case_when(str_detect(TRCODE,"^0101")~ "Sleeping",
                              str_detect(TRCODE,"^0102") ~ "Personal Activities",
                              str_detect(TRCODE,"^0103") ~ "Personal Activities",
                              str_detect(TRCODE,"^0104") ~ "Personal Activities",
                              str_detect(TRCODE,"^0105") ~ "Personal Activities",
                              str_detect(TRCODE,"^0199") ~ "Personal Activities",
                              str_detect(TRCODE,"^02") ~ "Household Activities",
                              str_detect(TRCODE,"^03") ~ "Helping HH Members",
                              str_detect(TRCODE,"^05") ~ "Work",
                              str_detect(TRCODE,"^06") ~ "School",
                              str_detect(TRCODE,"^07") ~ "Shopping",
                              str_detect(TRCODE,"^11") ~ "Eating and Drinking",
                              str_detect(TRCODE,"^12") ~ "Socializing and Leisure",
                              str_detect(TRCODE,"^13") ~ "Sports, Exercise, Recreation",
                              str_detect(TRCODE,"^14") ~ "Religious/Spiritual Activities",
                              str_detect(TRCODE,"^16") ~ "Telephone Calls",
                              str_detect(TRCODE,"^18") ~ "Traveling/Commuting",
                              TRUE ~ "Other"))%>%
  mutate(TUACTDUR24 = as.numeric(TUACTDUR24),
         startact = hms::as_hms(TUSTARTTIM),
         stopact = startact,
         stopact = stopact + hms::hms(minutes = TUACTDUR24),
         stopact = hms::hms(seconds = as.numeric(stopact)),
         stopact = ifelse(stopact > hms::as_hms("24:00:00"),
                          as.numeric(stopact) - (24*60*60),
                          as.numeric(stopact)),
         stopact = hms::hms(seconds = stopact),
         activity = as.factor(activity))%>%
  group_by(TUCASEID, activity, group_id = consecutive_id(activity)) %>%
  summarise(TUACTDUR24 = sum(TUACTDUR24),
            startact = min(startact),
            stopact = max(stopact),
            TUACTIVITY_N = mean(as.numeric(TUACTIVITY_N)), .groups = "drop") %>%
  arrange(group_id) %>%
  select(-group_id)

## ---- Additional Cleaning / Filtering ----
## Filter to not napping full time and part time students
no_nap_ftids <- diary_ftc %>%
  dplyr::filter(activity == "Sleeping")%>%
  group_by(TUCASEID)%>%
  arrange(startact)%>%
  mutate(instance = row_number())%>%
  slice_max(instance)%>%
  dplyr::filter(instance==2)%>%
  dplyr::select(TUCASEID)

no_nap_ptids <- diary_ptc %>%
  dplyr::filter(activity == "Sleeping")%>%
  group_by(TUCASEID)%>%
  arrange(startact)%>%
  mutate(instance = row_number())%>%
  slice_max(instance)%>%
  dplyr::filter(instance==2)%>%
  dplyr::select(TUCASEID)

diary_ftc <- diary_ftc %>%
  dplyr::filter(TUCASEID %in% no_nap_ftids$TUCASEID)
diary_ptc <- diary_ptc %>%
  dplyr::filter(TUCASEID %in% no_nap_ptids$TUCASEID)

## Filter to folks who had 2-3 meals 
meal_count_fc <- diary_ftc %>% 
  dplyr::filter(activity == "Eating and Drinking") %>% 
  group_by(TUCASEID) %>% 
  summarize(neach=n())%>%
  group_by(neach) %>% 
  summarize(n=n()) %>%
  mutate(p = round(n/sum(n),2))
# 68 students had 3 meals 139 had 2, 

meal_count_pc <- diary_ptc %>% 
  dplyr::filter(activity == "Eating and Drinking") %>% 
  group_by(TUCASEID) %>% 
  summarize(neach=n())%>%
  group_by(neach) %>% 
  summarize(n=n()) %>%
  mutate(p = round(n/sum(n),2))
# 32 students had 2 meals, 21 had 3 meals

id_meals_ft <- diary_ftc %>%
  dplyr::filter(activity == "Eating and Drinking")%>%
  arrange(TUACTIVITY_N)%>%
  group_by(TUCASEID)%>%
  mutate(mealno = row_number())%>%
  slice_max(mealno)%>%
  dplyr::filter(mealno %in% c(2,3))%>%
  dplyr::pull(TUCASEID)

id_meals_pt <- diary_ptc %>%
  dplyr::filter(activity == "Eating and Drinking")%>%
  arrange(TUACTIVITY_N)%>%
  group_by(TUCASEID)%>%
  mutate(mealno = row_number())%>%
  slice_max(mealno)%>%
  dplyr::filter(mealno %in% c(2,3))%>%
  dplyr::pull(TUCASEID)

diary_ftc <- diary_ftc %>%
  dplyr::filter(TUCASEID %in% id_meals_ft) 
# 207 Students
diary_ptc <- diary_ptc %>%
  dplyr::filter(TUCASEID %in% id_meals_pt)
# 53 students

write_rds(diary_ftc,here("data/proc_data","01_FTS.rds"))
write_rds(diary_ptc,here("data/proc_data","01_PTS.rds"))

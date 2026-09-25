## Simulations
# August 2026
# Sophia Arabadjis

library(here)
library(tidyverse)
library(simsurv)
library(flexsurv)

source(here("code","00_Source.R"))

dat <- readRDS(file = here("Output","03_TPM.rds"))

params <- read_csv(here("data/proc_data","02_SleepMealParams.csv"))

actlist <- read_csv(here("data/proc_data","03_ActivityList.csv"))
actlistsmall <- paste(actlist$activity[-c(1:3,11)])

# trial <- sleep_sim.f(days=2,num_hh=10,ind_hh=2,group="Full-Time",params=params)
# trial2 <- meal_sim.f(days=2,num_hh=10,ind_hh=2,group="Full-Time",iact = trial)
# trial3 <- crisk.f(group="Full-Time",iact2 = rbind(trial,trial2),num_hh = 10, ind_hh = 2,days = 2,
#                   actlist = actlistsmall)

FT_sim <- big_sim.f(days = 3, num_hh = 200, ind_hh=5, group="Full-Time",actlist=actlistsmall,params=params)
PT_sim <- big_sim.f(days = 3, num_hh = 200, ind_hh=5, group="Part-Time",actlist=actlistsmall,params=params)

write_rds(FT_sim,here("Output","04_FuTiSim.rds"))
write_rds(PT_sim,here("Output","04_PaTiSim.rds"))

# Compute details

profile_runsleep <- compute_details.f(
    func     = sleep_sim.f, 
    group    = "Full-Time", 
    num_hh   = 200, 
    ind_hh   = 5, 
    days     = 3,
    params  = params
)

profile_runmeal <- compute_details.f(
    func     = meal_sim.f, 
    group    = "Full-Time", 
    num_hh   = 200, 
    ind_hh   = 5, 
    days     = 3,
    iact = profile_runsleep$output
)

profile_run <- compute_details.f(
    func     = crisk.f, 
    group    = "Full-Time", 
    num_hh   = 200, 
    ind_hh   = 5, 
    days     = 3, 
    iact2 = profile_runmeal$output,
    actlist  = colnames(dat$FTP)
    )

profile_big <-compute_details.f(
    func = big_sim.f,
    group = "Full-Time",
    num_hh = 200,
    ind_hh = 5,
    days = 3,
    actlist=actlistsmall,
    params=params
)

compute_res <- data.frame(full = unlist(profile_big$metrics),
                          sleep = unlist(profile_runsleep$metrics),
                          meal = unlist(profile_runmeal$metrics),
                          crisk = unlist(profile_run$metrics))%>%
    rownames_to_column(var = "Function")

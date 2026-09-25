# HASP
Hierarchical Activity Simulation Procedure

A procedure for modeling and simulating daily time-use schedules of full-time and
part-time college students, an example built on American Time Use Survey (ATUS) data.
The project fits distributions to sleep and meal timing, estimates activity transition
probabilities, and generates synthetic, minute-resolution activity diaries via a
competing-risks simulation framework. A validation of the simulated schedules
against the observed ATUS data is provided to match the paper by Arabadjis and Sweeney,
2026.

Arabadjis and Sweeney. ``Hierarchical Activity Simulation Procedure (HASP): 
A Time-to-event-based method for generating realistic activity patterns." 2026. 

## Overview

The procedure has three stages:

1. **Data processing** — Read raw ATUS `.dat` files, filter to full-time and
   part-time students, recode detailed activity codes into broader categories, and
   clean/filter diaries (dropping nappers, keeping students with 2–3 recorded meals).
2. **Parameter estimation** — Fit statistical distributions to sleep duration and
   inter-meal timing (including a custom zero-inflated exponential/normal family for
   "time to first meal," since many students eat immediately upon waking), and
   compute activity-transition probability matrices (TPMs), including period-specific
   TPMs for the gaps between sleep and meals.
3. **Simulation & validation** — Simulate correlated household sleep schedules
   (clustered by household), simulate meal timing from the fitted distributions, fill
   remaining time with a competing-risks draw over transition probabilities, and
   compare simulated activity profiles/survival curves against the empirical ATUS data.

## Requirements

R with the following packages:

```r
install.packages(c(
  "here", "tidyverse", "survival", "survminer", "flexsurv",
  "gamlss.dist", "fitdistrplus", "simsurv", "patchwork", "hms", "stringr"
))
```

## Project Structure

```
code/
  00_Source.R              Shared helper functions (see below) — sourced by every downstream script
  01_ATUSprocessing.R       Cleans raw ATUS files into analysis-ready student diaries (01_FTS.rds, 01_PTS.rds)
  02_ParametersStatDist.R   Fits sleep/meal timing distributions; writes 02_SleepMealParams.csv, 02_ActDist.csv
  03_TPM.R                  Builds transition probability matrices and period-based TPMs; writes 03_TPM.rds
  04_Simulations.R          Runs the full household simulation (big_sim.f); profiles runtime/memory
  05_SimResults.R           Post-processes simulated diaries into minute-by-minute activity profiles and figures
  06_ATUScheck.R            Builds the same minute-by-minute profiles from observed ATUS data for comparison
data/
  atusact-2024/, atusrost-2024/, atusresp-2024/   Raw ATUS activity/roster/respondent files (not included)
  proc_data/                Cleaned diaries and fitted parameters written by the pipeline
Output/
  Simulation results, transition matrices, and generated figures
```

Paths are resolved with the `here` package, so scripts should be run from the
project root (the folder containing `code/` and `data/`).

## Procedure / Run Order

```r
source("code/01_ATUSprocessing.R")     # -> data/proc_data/01_FTS.rds, 01_PTS.rds
source("code/02_ParametersStatDist.R") # -> data/proc_data/02_SleepMealParams.csv, 02_ActDist.csv
source("code/03_TPM.R")                # -> Output/03_StationaryDistributions.csv, Output/03_TPM.rds
source("code/04_Simulations.R")        # -> Output/04_FuTiSim.rds, 04_PaTiSim.rds
source("code/05_SimResults.R")         # simulated activity profiles, survival curves, figures
source("code/06_ATUScheck.R")          # observed (ATUS) activity profiles, for comparison to 05
```

`00_Source.R` is sourced automatically inside scripts `02`–`06` and does not need
to be run on its own. 

Raw ATUS `.dat` files (activity, roster, and respondent files) are expected under
`data/atusact-2024/`, `data/atusrost-2024/`, and `data/atusresp-2024/` — these are
external ATUS extracts and are not included in the repo.

## Key Functions (`00_Source.R`)

**Zero-inflated distributions** — custom `d`/`p`/`q`/`r` functions (in the standard
R distribution-function style) used to model "time to first meal," where a
meaningful fraction of students eat immediately after waking:
- `dZIE`, `pZIE`, `qZIE`, `rZIE` — zero-inflated exponential
- `dZIN`, `pZIN`, `qZIN`, `rZIN` — zero-inflated normal
- `dTF_safe`, `pTF_safe`, `qTF_safe`, `rTF_safe` — parameter-validated wrappers
  around `gamlss.dist`'s t-distribution family

**Transition matrices**
- `tm.f()` — builds one individual's activity transition-count matrix
- `tm_tab.f()` — aggregates transition matrices across individuals into a
  row-normalized transition probability matrix, iterated to a stationary product
- `tm_Period_tab.f()` — transition probabilities restricted to a period (e.g.
  between sleep and first meal)
- `period_diaries.f()` — splits each diary into the three periods used for
  period-based transition probabilities (pre-meal-1, between meals, post-last-meal)

**Sampling**
- `SVP.f()`, `samp_rnd.f()`, `propz()`, `activity.f()` — supporting sampling
  utilities (scaled probability vectors, fractional-count rounding, proportion
  confidence intervals)

**Simulation** (the core generative model, run per household group)
- `sleep_sim.f()` — simulates household-clustered sleep episodes using
  `simsurv`, fit to a Weibull (part-time) or Gaussian hazard (full-time) model,
  but modifiable for other distributions.
- `meal_sim.f()` — simulates meal timing relative to sleep, using the fitted
  zero-inflated/gamma/normal meal-gap distributions (again, modifiable for
  other distributions.)
- `crisk.f()` — fills remaining gaps between sleep and meals using a
  competing-risks draw over period-based transition probabilities
- `big_sim.f()` — orchestrates `sleep_sim.f()` → `meal_sim.f()` → `crisk.f()`
  into one complete simulated diary set
- `compute_details.f()` — profiles a function's runtime and peak memory use for
  a single call (used to benchmark the simulation steps)

**Cleaning / visuals**
- `daysect.f()` — splits simulated activity spans that cross midnight into
  same-day segments, for plotting/analysis
- `surv_plots.f()` — builds survival-analysis-ready data for a given activity
  and occurrence number (e.g. "2nd meal of the day"), handling censoring for
  individuals who don't reach that occurrence

## Notes

- `%notin%` is defined as the negation of `%in%` and used throughout.
- Activity codes are recoded from ATUS `TRCODE` prefixes into categories such as
  Sleeping, Personal Activities, Household Activities, Work, School, Eating and
  Drinking, Socializing/Leisure, and Other; a handful of categories (e.g. legal
  services, government/civic service, volunteering) are folded into "Other."
- All simulation times are in minutes, with day boundaries at multiples of 1440.

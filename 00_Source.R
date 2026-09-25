## Source Code
# August 2026
# Sophia Arabadjis

# ---- Useful Functions ----

##############################
###### General Functions #####
##############################
`%notin%` <- Negate(`%in%`)


# ---- Zero Inflated Exponential Suite ----

# 1. Probability Density Function (PDF)
dZIE <- function(x, rate = 1, pi = 0.1, log = FALSE) {
  if (rate <= 0) stop("Rate must be positive.")
  if (pi < 0 || pi > 1) stop("Pi must be between 0 and 1.")
  
  density <- rep(0, length(x))
  
  # Structural zeros get the discrete probability mass 'pi'
  density[x == 0] <- pi
  
  # Positive values follow the scaled continuous exponential density
  pos_idx <- x > 0
  density[pos_idx] <- (1 - pi) * dexp(x[pos_idx], rate = rate)
  
  if (log) return(log(density))
  return(density)
}

# 2. Cumulative Distribution Function (CDF)
pZIE <- function(q, rate = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE) {
  if (rate <= 0) stop("Rate must be positive.")
  if (pi < 0 || pi > 1) stop("Pi must be between 0 and 1.")
  
  p <- rep(0, length(q))
  
  # Values strictly below zero have no probability mass
  p[q < 0] <- 0
  
  # At or above zero, add structural zeros 'pi' to the scaled exponential CDF
  ge_idx <- q >= 0
  p[ge_idx] <- pi + (1 - pi) * pexp(q[ge_idx], rate = rate)
  
  if (!lower.tail) p <- 1 - p
  if (log.p) return(log(p))
  return(p)
}

# 3. Quantile Function (Inverse CDF)
qZIE <- function(p, rate = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE) {
  if (rate <= 0) stop("Rate must be positive.")
  if (pi < 0 || pi > 1) stop("Pi must be between 0 and 1.")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("Probabilities must be between 0 and 1.")
  
  q <- rep(NA, length(p))
  
  # Any probability falling within the structural zero mass maps directly to 0
  q[p <= pi] <- 0
  
  # Probabilities above 'pi' map to the shifted/scaled exponential quantiles
  gt_idx <- p > pi
  q[gt_idx] <- qexp((p[gt_idx] - pi) / (1 - pi), rate = rate)
  
  return(q)
}

# 4. Random Generation Function
rZIE <- function(n, rate = 1, pi = 0.1) {
  if (rate <= 0) stop("Rate must be positive.")
  if (pi < 0 || pi > 1) stop("Pi must be between 0 and 1.")
  
  # Determine which indices are structural zeros
  is_zero <- rbinom(n, size = 1, prob = pi)
  exp_values <- rexp(n, rate = rate)
  
  return(ifelse(is_zero == 1, 0, exp_values))
}

# ---- Zero Inflated Normal Suite ----
# 1. Density Function
dZIN <- function(x, mean = 0, sd = 1, pi = 0.1, log = FALSE) {
  if (sd <= 0) stop("sd must be positive.")
  if (pi < 0 || pi > 1) stop("pi must be between 0 and 1.")
  density <- rep(0, length(x))
  density[x == 0] <- pi
  non_zero_idx <- x != 0
  density[non_zero_idx] <- (1 - pi) * dnorm(x[non_zero_idx], mean = mean, sd = sd)
  if (log) return(log(density))
  return(density)
}

# 2. Cumulative Distribution Function
pZIN <- function(q, mean = 0, sd = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE) {
  if (sd <= 0) stop("sd must be positive.")
  if (pi < 0 || pi > 1) stop("pi must be between 0 and 1.")
  p <- rep(0, length(q))
  
  # Below zero
  neg_idx <- q < 0
  p[neg_idx] <- (1 - pi) * pnorm(q[neg_idx], mean = mean, sd = sd)
  
  # At or above zero (Fixed brace from prior version)
  ge_idx <- q >= 0
  p[ge_idx] <- pi + (1 - pi) * pnorm(q[ge_idx], mean = mean, sd = sd)
  
  if (!lower.tail) p <- 1 - p
  if (log.p) return(log(p))
  return(p)
}

# 3. Quantile Function
qZIN <- function(p, mean = 0, sd = 1, pi = 0.1, lower.tail = TRUE, log.p = FALSE) {
  if (sd <= 0) stop("sd must be positive.")
  if (pi < 0 || pi > 1) stop("pi must be between 0 and 1.")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("Probabilities must be between 0 and 1.")
  
  q <- rep(NA, length(p))
  p_zero_norm <- pnorm(0, mean = mean, sd = sd)
  cdf_at_zero_minus <- (1 - pi) * p_zero_norm
  cdf_at_zero_plus  <- pi + (1 - pi) * p_zero_norm
  
  # Assign regions
  q[p < cdf_at_zero_minus] <- qnorm(p[p < cdf_at_zero_minus] / (1 - pi), mean = mean, sd = sd)
  q[p >= cdf_at_zero_minus & p <= cdf_at_zero_plus] <- 0
  q[p > cdf_at_zero_plus] <- qnorm((p[p > cdf_at_zero_plus] - pi) / (1 - pi), mean = mean, sd = sd)
  return(q)
}

# 4. Random Generation Function
rZIN <- function(n, mean = 0, sd = 1, pi = 0.1) {
  if (sd <= 0) stop("sd must be positive.")
  if (pi < 0 || pi > 1) stop("pi must be between 0 and 1.")
  samples <- rnorm(n, mean = mean, sd = sd)
  is_zero <- as.logical(rbinom(n, size = 1, prob = pi))
  samples[is_zero] <- 0
  return(samples)
}

# Functions for T-Distribution
# 1. Create a safe wrapper for the density function (dTF)
dTF_safe <- function(x, mu = 0, sigma = 1, nu = 1, log = FALSE) {
  # If parameters are invalid, return NaN vectors matching the input length
  if (sigma <= 0 || nu <= 0) {
    return(rep(NaN, length(x)))
  }
  # Otherwise, use the standard gamlss function
  gamlss.dist::dTF(x, mu = mu, sigma = sigma, nu = nu, log = log)
}

# 2. Create a safe wrapper for the cumulative distribution function (pTF)
pTF_safe <- function(q, mu = 0, sigma = 1, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  if (sigma <= 0 || nu <= 0) {
    return(rep(NaN, length(q)))
  }
  gamlss.dist::pTF(q, mu = mu, sigma = sigma, nu = nu, lower.tail = lower.tail, log.p = log.p)
}

# 3. Quantile Function Wrapper
qTF_safe <- function(p, mu = 0, sigma = 1, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  # If log.p is TRUE, p can be negative. If FALSE, p must be between 0 and 1.
  if (sigma <= 0 || nu <= 0) return(rep(NaN, length(p)))
  if (!log.p && any(p < 0 | p > 1, na.rm = TRUE)) return(rep(NaN, length(p)))
  
  gamlss.dist::qTF(p, mu = mu, sigma = sigma, nu = nu, lower.tail = lower.tail, log.p = log.p)
}

# 4. Random Generation Function Wrapper
rTF_safe <- function(n, mu = 0, sigma = 1, nu = 1) {
  if (sigma <= 0 || nu <= 0) return(rep(NaN, n))
  gamlss.dist::rTF(n, mu = mu, sigma = sigma, nu = nu)
}

##############################
### Functions for Visuals ###
##############################

surv_plots.f <- function(v,n){
  t <- rbind(diary_ftc %>%
               mutate(Student = "Full-Time"),
             diary_ptc %>%
               mutate(Student = "Part-Time"))
  m <- t %>%
    dplyr::filter(activity == v)%>%
    group_by(TUCASEID)%>%
    mutate(number = row_number())%>%
    dplyr::filter(number == n)
  if(length(which(t$TUCASEID %notin% m$TUCASEID))==0){
    temp <- t %>%
      dplyr::filter(activity == v)%>%
      group_by(TUCASEID)%>%
      mutate(number = row_number())%>%
      dplyr::filter(number == n)%>%
      mutate(min_elapse = hms::as_hms(startact - hms::as_hms("00:00:00")),
             dur_elapse = as.numeric(min_elapse)/(60*60),
             status=1)%>%
      ungroup()
  }
  if(length(which(t$TUCASEID %notin% m$TUCASEID))>0){
    censored <- t %>%
      dplyr::select(TUCASEID,Student)%>%
      distinct()%>%
      mutate(status = ifelse(TUCASEID %notin% m$TUCASEID,0,1),
             activity = v,
             number = n,
             dur_elapse = 20)%>%
      dplyr::filter(status == 0)
    
    temp <- t %>%
      dplyr::filter(activity == v)%>%
      group_by(TUCASEID)%>%
      mutate(number = row_number())%>%
      dplyr::filter(number == n)%>%
      mutate(min_elapse = hms::as_hms(startact - hms::as_hms("00:00:00")),
             dur_elapse = as.numeric(min_elapse)/(60*60),
             status = 1)%>%
      ungroup()%>%
      dplyr::select(c(TUCASEID,status,dur_elapse,number,Student,activity))
    
    temp <- rbind(temp,censored)
  }
  
  fit <- survfit(Surv(dur_elapse,event = status)~Student,data = temp) 
  
  fig <- ggsurvplot(fit,
                    palette = c('cadetblue','orange'),
                    conf.int = T,
                    pval=F,
                    title = paste(v,": Event ",n,sep=""),
                    surv.median.line = "hv",
                    data = temp)
  
  return(list(fit,fig))
  
}

# ---- Simple Sampling -----

# Define activity replace function
activity.f <- function(v, a) {
  v <- ifelse(v == 1 | v == "1", a, v)
}

# Define scaled probability vector
# (because probabilities increase for each binom after x category is removed)
SVP.f <- function(probvector) {
  SVP <- NA
  SVP[1] <- probvector[1]
  for (k in 2:length(probvector)) {
    SVP[k] <- probvector[k] / (1 - sum(probvector[1:k - 1]))
  }
  return(SVP)
}

propz <- function(value,samplesize,alpha=.05){
  qnorm(1-alpha/2)*sqrt((value*(1-value))/samplesize)
}

# Function for simulating from a Multinomial distribution using binomials
## Random Sampler, No recall Per Person Aggregation
MNsim_sampler_P.f <- function(nsize, tsize, tvec, probvector, activityvector, output = "results",u = 'minutes') {
  if (u == 'minutes') {
    mins <- 1440
  } else if (is.numeric(u)) {
    mins <- u
  } else {
    stop("Argument 'u' must be 'minutes' or a specific numeric value.")
  }
  mat <- matrix(data = 0,
                ncol = nsize,
                nrow = mins) # Change to minute time chunks
  SVP <-round(SVP.f(probvector = probvector),5)
  # Simulation
  for (j in 1:dim(mat)[2]) {
    for (i in 1:dim(mat)[1]) {
      for (k in 1:length(SVP)) {
        mat[i, j] <- ifelse(mat[i, j] == 0, rbinom(1, 1, SVP[k]), mat[i, j])
        mat[, j] <- activity.f(mat[, j], activityvector[k])
      }
    }
  }
  # Random sampler, 3 
  mat_sample_R <- mat[c(sample(18,replace=F,size=tsize)),]
  # Fixed sampler
  mat_sample_F <- mat[tvec,]
  # Output
  out_R <- data.frame(mat_sample_R) %>%
    pivot_longer(cols = starts_with("X"),
                 names_to = "SimPerson",
                 values_to = "activity") %>%
    group_by(SimPerson, activity) %>%
    summarise(obs_count = n()) %>%
    ungroup() %>%
    complete(nesting(SimPerson),
             nesting(activity),
             fill = list(obs_count = 0))%>%
    group_by(SimPerson)%>%
    mutate(obs_prop = obs_count / sum(obs_count)) %>%
    ungroup() %>%
    group_by(activity) %>%
    summarize(mean_obs_prop = mean(obs_prop)) %>%
    arrange(activity)%>%
    mutate(activity = factor(activity,levels=c(1:14),labels=actlist))%>%
    complete(activity,fill = list(mean_obs_prop = 0))
  out_F <- data.frame(mat_sample_F) %>%
    pivot_longer(cols = starts_with("X"),
                 names_to = "SimPerson",
                 values_to = "activity") %>%
    group_by(SimPerson, activity) %>%
    summarise(obs_count = n()) %>%
    ungroup() %>%
    complete(nesting(SimPerson),
             nesting(activity),
             fill = list(obs_count = 0))%>%
    group_by(SimPerson)%>%
    mutate(obs_prop = obs_count / sum(obs_count)) %>%
    ungroup() %>%
    group_by(activity) %>%
    summarize(mean_obs_prop = mean(obs_prop)) %>%
    arrange(activity)%>%
    mutate(activity = factor(activity,levels=c(1:14),labels=actlist))%>%
    complete(activity,fill = list(mean_obs_prop = 0))
  # out <- data.frame(fixed = out_F$mean_obs_prop,
  #                   random = out_R$mean_obs_prop,
  #                   activity = out_F$activity)
  if(output == "samplematrix"){return(mat)}
  #if(output == "results"){return(out)}
  if(output == "verbose"){return(list(out,mat))}
}


## ---- Transition Matrices ----

tm.f <- function(df, id_list, id,acts){
  key <- data.frame(label = paste(actlist),level=1:length(actlist))
  key <- key %>%
    dplyr::filter(label %in% acts)
  # Select individuals
  pd <- df[df$TUCASEID == id_list[id],]
  # Tabulate diagonals
  pm <- matrix(data=0,nrow=length(acts),ncol=length(acts))
  diags <- pd  %>%
    group_by(activity) %>% 
    summarize(mins=sum(TUACTDUR24),
              n=n()) %>%
    mutate(n = ifelse(is.na(n),0,n),
           mins = mins-n)%>%
    complete(activity,fill=list(mins=0))
  diag(pm) <- diags$mins[which(diags$activity %in% key$label)]
  colnames(pm) <- acts
  rownames(pm) <- acts
  
  # Tabulate transitions
  pdlong <- data_frame(from = pd$activity[-dim(pd)[1]],
                       to= factor(c(paste(pd$activity[-1])),
                                  levels = levels(pd$activity)),
                       t = 1)%>%
    complete(from,to,fill=list(t=0))
  pdlong <- pdlong %>%
    summarize(t=sum(t),.by=c(from,to))%>%
    arrange(from,to)%>%
    pivot_wider(names_from = to,values_from = t)%>%
    dplyr::select(sort(names(.)))%>%
    dplyr::select(-from)
  # Sum to create full TM
  if(length(acts) < length(actlist)){
    pdlong <- pdlong[c(key$level),c(key$level)]
  }
  pm <- pm+pdlong
  # pm2 <- as.matrix(pm/rowSums(pm))
  # pm2[is.nan(pm2)] <- 0
  return(pm)
}

period_diaries.f <- function(d){
  actno_meal1 <- d %>%
    dplyr::filter(activity=="Eating and Drinking 1",.by = TUCASEID)%>%
    mutate(locs = TUACTIVITY_N)%>%
    dplyr::select(TUCASEID,locs)
  actno_maxmeal <- d %>%
    dplyr::filter(act == "Eating and Drinking" & 
                    countnum == max(countnum),.by=TUCASEID)%>%
    mutate(locs2 = TUACTIVITY_N)%>%
    dplyr::select(TUCASEID,locs2)
  Period1 <- d %>%
    left_join(actno_meal1)%>%
    dplyr::filter(TUACTIVITY_N < locs)%>%
    dplyr::filter(act != "Sleeping")
    Period2 <- d %>%
      left_join(actno_meal1)%>%
      dplyr::filter(TUACTIVITY_N > locs)%>%
      left_join(actno_maxmeal)%>%
      dplyr::filter(TUACTIVITY_N < locs2)
    Period3 <- d %>%
      left_join(actno_maxmeal)%>%
      dplyr::filter(TUACTIVITY_N > locs2)%>%
      dplyr::filter(act != "Sleeping")
    return(list(Period1 = Period1,Period2 = Period2,Period3 = Period3))
}

# ---- Full Transition Matrix Function ----

tm_tab.f <- function(df,id_list,acts = actlist,id,verbose=FALSE){
  # Matrix
  tm <- array(unlist(lapply(1:length(id_list), 
                            tm.f, 
                            df= df,
                            id_list=id_list,
                            acts = acts)), 
              dim = c(length(acts),length(acts),length(id_list)))
  tmP <- apply(tm, c(1,2), FUN = sum)
  tmP <- as.matrix(tmP/rowSums(tmP))
  colnames(tmP)<-acts
  rownames(tmP)<-acts
  
  # TMP
  TmP <- tmP1 <- tmP
  dft <- .01
  aft <- .01
  
  while(dft > 0.0001){
    dft <- mean(abs(as.vector(TmP-(TmP %*% tmP1))))
    aft <- c(aft,dft)
    TmP <- TmP %*% tmP1
  }
  
  if(verbose==FALSE){return(list(Raw = tmP,Product = TmP))}
  if(verbose==TRUE){return(list(Raw=tmP,Product=TmP,Trace = aft))}
}

tm_Period_tab.f <- function(df,id_list,acts = actlist,id,Period){
  # Matrix
  tm <- array(unlist(lapply(1:length(id_list), 
                            tm.f, 
                            df= df,
                            id_list=id_list,
                            acts = acts),), 
              dim = c(length(acts),length(acts),length(id_list)))
  tmP <- apply(tm, c(1,2), FUN = sum)
  rownames(tmP) <- acts
  colnames(tmP) <- acts
  tmCS <- colSums(tmP)  
  tmCS <- as.matrix(tmCS/sum(tmCS))
  rownames(tmCS)<-acts
  colnames(tmCS)<-paste("Period",Period,sep="_")
  tmCS <- tmCS[order(tmCS[,1],decreasing = T),]
  tmCS <- t(tmCS)
  
  return(tmCS)
}

## ---- Monte Carlo Sampling ----
# Function for randomly distributing fractional counts
samp_rnd.f<-function(d){ 
  tot<-sum(d)
  w_<-which(d>0)
  v_<-sample(w_,length(w_)-1)
  d[v_]<-round(d[v_])
  d[setdiff(w_,v_)]<-tot-sum(d[v_])
  return(d)
  }

# Per person approach, no Multinomial For now

MN_corr_sampler.f <- function(mm, av, nsize, tsize, tvec, reps, output = "samplematrix"){
  
  radix <- 10000 
  rep <- reps*10 
  
  # Initialize matrix of individual activities
  iact <- matrix(NA, nrow = radix, ncol = rep) 
  iact[, 1] <- av[1] # Start everyone explicitly in the first activity state
  
  # SAFE SAMPLING FUNCTION: Prevents R from breaking when a pool has only 1 person
  safe_sample <- function(x, size, ...) {
    if (length(x) <= 1) {
      return(rep(x, size))
    } else {
      return(sample(x, size, ...))
    }
  }
  
  # Simulation Loops
  for(i in 2:rep){
    # Create a temporary container for time step i so updates don't bleed into each other
    next_step_vector <- iact[, i - 1]
    
    for(j in 1:(length(av))){
      
      # Find individuals who are in activity j at the end of the previous time step
      w <- which(iact[, i - 1] == av[j])
      
      if(length(w) > 0){
        
        # Calculate exactly how many people should transition based on mm rows summing to 1
        expected_counts <- length(w) * mm[j, ]
        
        # Call your exact rounding function
        d <- samp_rnd.f(expected_counts) 
        
        for(k in 1:length(d)){
          # SAFETY FIX: Cap d[k] to ensure it never asks for more elements than what's left in w
          # Also ensures negative allocations generated by samp_rnd.f are ignored
          allocation_size <- max(0, min(d[k], length(w)))
          
          if(allocation_size > 0){
            # Use safe_sample instead of sample()
            v <- safe_sample(w, allocation_size)
            
            # Update our temporary vector instead of mutating the tracking matrix live
            next_step_vector[v] <- av[k]
            w <- setdiff(w, v)
          }
          
          # Your kluge fix for residual allocations
          if(k == length(d) & length(w) > 0){
            next_step_vector[w] <- av[k]
          } 
        }
      }
    }
    # Commit all transitions simultaneously for time step i
    iact[, i] <- next_step_vector
  }
  
  # Corrected row sampling mechanism 
  sampled_rows <- sample(1:nrow(iact), size = nsize, replace = FALSE)
  mat <- t(iact[sampled_rows, c((9*reps):(10*reps))])
  
  return(mat)
}


MN_corr_sampler2.f <- function(mm, av, nsize, tsize, tvec, output = "samplematrix"){
  # mm = markhov matrix / transition matrix 
  # av = activity vector / names of activities corresponding to Markov Matrix
  # nsize = # of individuals
  # tsize = # of samples for random sampler
  # tvec = # and location of fixed sampler (indexed 1:48)
  # iact = individual activity matrix (sample matrix)
  
  radix<-10000 # initial population size; all in activity "A" which is sleeping
  rep<-300 # discrete replications
  P<-c(radix,rep(0,14)) # initial distribution of activities
  iact<-matrix(NA,radix,rep) # initialize matrix of individual activities; individuals on rows, replications across columns
  
  iact[,1]<-av[1]
  for(i in 2:rep){
    for(j in 1:(length(av))){
      w<-which(iact[,i-1]==av[j])
      if(length(w)>0){
        d<-samp_rnd.f(P[j]%*%mm[j,]) # could incorporate multinomial sampling in this line; currently just element x vector of tran matrix
        for(k in 1:length(d)){
          if(d[k]>0){
           # print(paste(i,j,k,d[k]))
            v<-sample(w,d[k])
            iact[v,i]<-av[k]
            w<-setdiff(w,v)}
          if(k==length(d) & length(w)>0){iact[w,i]<-av[k]} # kluge fix for an error; residual w results in NA, allocate to activity[last]
        }
      }}
    P<-sapply(av,function(x) sum(iact[,i]==x))
  }
  mat <- t(iact[sample(dim(iact)[1],size = nsize,replace=F),c(253:300)])
  
  # Random sampler 
  mat_sample_R <- mat[c(sample(48,replace=F,size=tsize)),]
  # Fixed sampler
  mat_sample_F <- mat[tvec,]
  # Output
  out_R <- data.frame(mat_sample_R) %>%
    pivot_longer(cols = starts_with("X"),
                 names_to = "SimPerson",
                 values_to = "activity") %>%
    group_by(SimPerson, activity) %>%
    summarise(obs_count = n()) %>%
    ungroup() %>%
    complete(nesting(SimPerson),
             nesting(activity),
             fill = list(obs_count = 0))%>%
    group_by(SimPerson)%>%
    mutate(obs_prop = obs_count / sum(obs_count)) %>%
    ungroup() %>%
    group_by(activity) %>%
    summarize(mean_obs_prop = mean(obs_prop)) %>%
    arrange(activity) # Alphabetizes Activities
  out_F <- data.frame(mat_sample_F) %>%
    pivot_longer(cols = starts_with("X"),
                 names_to = "SimPerson",
                 values_to = "activity") %>%
    group_by(SimPerson, activity) %>%
    summarise(obs_count = n()) %>%
    ungroup() %>%
    complete(nesting(SimPerson),
             nesting(activity),
             fill = list(obs_count = 0))%>%
    group_by(SimPerson)%>%
    mutate(obs_prop = obs_count / sum(obs_count)) %>%
    ungroup() %>%
    group_by(activity) %>%
    summarize(mean_obs_prop = mean(obs_prop)) %>%
    arrange(activity)
  out_T <-data.frame(mat) %>%
    pivot_longer(cols = starts_with("X"),
                 names_to = "SimPerson",
                 values_to = "activity") %>%
    group_by(SimPerson, activity) %>%
    summarise(obs_count = n()) %>%
    ungroup() %>%
    complete(nesting(SimPerson),
             nesting(activity),
             fill = list(obs_count = 0))%>%
    group_by(SimPerson)%>%
    mutate(obs_prop = obs_count / sum(obs_count)) %>%
    ungroup() %>%
    group_by(activity) %>%
    summarize(mean_obs_prop = mean(obs_prop)) %>%
    arrange(activity)
  # out <- data.frame(Fixed = out_F$mean_obs_prop,
  #                   Random = out_R$mean_obs_prop,
  #                   Unobserved = out_T$mean_obs_prop,
  #                   Activity = sort(out_F$activity))
  if (output == "samplematrix") {return(mat)}
  if (output == "sampleresults") {return(out)}
  #all <- list(mat,out)
  #if (output == "verbose") {return(all)}
}

# ----------------------
# # Sleep Function # # 
# ----------------------

# Purpose: generate sleep correlated sleep functions for individuals within a 
# household. Ultimately nested within a broader function.

sleep_sim.f <- function(days, num_hh, ind_hh, group, params){
  if (group == 'Part-Time') {
    group = "Part-Time"
  } else if (group == "Full-Time") {
    group = "Full-Time"
  } else {
    stop("Argument 'group' must be either Full-Time or Part-Time. Else modify function.")
  }
  if(class(days) == "character"){stop("Days must be numeric.")}
  if(class(num_hh) == "character"){stop("Number of households must be numeric.")}
  if(class(ind_hh)== "character"){stop("Number of individuals per household must be numeric.")}
  
  # num_hh = number of households
  # num_ind = number of individuals in each household (set size)
  # sample_size = number of total individuals
  
  num_hh <- num_hh
  num_ind <- ind_hh
  samplesize <- num_hh * ind_hh
  
  # Covariate Data
  # hh = household number 
  # sex = arbitrary, can modify (can also introduce additional covariates here)
  covdat_sleep <- data.frame(id = 1:samplesize, hh = rep(1:num_hh,each = num_ind),
                             sex = rep(c("M","F"), length.out=samplesize))
  
  # Group Parameters and Initiate Survival Times:
  if(group == "Part-Time"){
    g <-params$Params[params$Pop=="Part-Time" & params$Var == "Sleep" & params$Labs=="shape"]
    l <- 1/params$Params[params$Pop=="Part-Time" & params$Var == "Sleep" & params$Labs=="scale"]^g
      
    d <- simsurv(dist = "weibull", 
            x = covdat_sleep, 
            lambdas = l,
            gammas = g,
            maxt = NULL,
            interval = c(0,1440))
            #interval = c(1E-8, 100000))
  }
  if(group == "Full-Time"){
    mu = params$Params[params$Pop=="Full-Time" & params$Var == "Sleep" & params$Labs=="mean"]
    sigma = params$Params[params$Pop=="Full-Time" & params$Var == "Sleep" & params$Labs=="sd"]
    betas = c(mean=mu,sd=sigma)
    
    norm_haz.f <- function(x, mu,s){
      res <- dnorm(x,mean=mu,sd=s)/pnorm(x,mean=mu,sd=s,lower.tail = F)
    }
    haz_sleep.f <- function(t,x,betas){
      res <- norm_haz.f(x = t,mu=betas[1],s=betas[2])
    }
    d <- simsurv(betas = betas, 
                 x = covdat_sleep, 
                 hazard = haz_sleep.f, 
                 maxt = NULL,
                 interval = c(0, 1440))
  }
  d_T <- d$eventtime
  
  # Cluster by household
  tstart_0 <- sample((-180:60)+1440, size = num_hh, replace = T)
  
  # arbitrarily set sleep sample range to be pre-midnight
  tstart_0 <- rep(tstart_0,each = num_ind)
  tstart_0 <- tstart_0 + rnorm(samplesize, mean = 0, sd = 15)
  tstop_0 <- tstart_0+d_T
  
  # Data Frame (individual activity = iact)
  iact = data.frame(id = 1:samplesize, HH = covdat_sleep$hh, sex = covdat_sleep$sex,
                    activity = "Sleeping", tstart = tstart_0, tstop=tstop_0,
                    event = 0)
  
  # Simulate other days
  for (i in 1:days){
    if(group == "Part-Time"){
      g <-params$Params[params$Pop=="Part-Time" & params$Var == "Sleep" & params$Labs=="shape"]
      l <- 1/params$Params[params$Pop=="Part-Time" & params$Var == "Sleep" & params$Labs=="scale"]^g
      
      d <- simsurv(dist = "weibull", 
                   x = covdat_sleep, 
                   lambdas = l,
                   gammas = g,
                   maxt = NULL,
                   interval = c(0,1440))
    }
    if(group == "Full-Time"){
      mu = params$Params[params$Pop=="Full-Time" & params$Var == "Sleep" & params$Labs=="mean"]
      sigma = params$Params[params$Pop=="Full-Time" & params$Var == "Sleep" & params$Labs=="sd"]
      betas = c(mean=mu,sd=sigma)
      
      norm_haz.f <- function(x, mu,s){
        res <- dnorm(x,mean=mu,sd=s)/pnorm(x,mean=mu,sd=s,lower.tail = F)
      }
      haz_sleep.f <- function(t,x,betas){
        res <- norm_haz.f(x = t,mu=betas[1],s=betas[2])
      }
      d <- simsurv(betas = betas, 
                   x = covdat_sleep, 
                   hazard = haz_sleep.f, 
                   maxt = NULL,
                   interval = c(0, 1440))
    }
    # Extract survival times
    d_T <- d$eventtime
    # Generate "start" times for sleep (Clustered)
    tstart <- sample((-180:60)+1440*(i+1), size = num_hh, replace=T)
    tstart <- rep(tstart,each = num_ind)
    # Add noise
    tstart <- tstart + rnorm(samplesize, mean = 0, sd = 15)
    #Stop Time
    tstop <- tstart+d_T
    temp <- data.frame(id = 1:samplesize, HH = covdat_sleep$hh, sex = covdat_sleep$sex, 
                       activity = "Sleeping", tstart, tstop, event=i)
    iact <- rbind(iact,temp)
  }
  return(iact)
}

# ----------------------
# # Meals Function # #
# ----------------------

# Purpose: generate meal times for individuals within a 
# household. Not currently clustered for common meals; 
# could introduce clustering. Ultimately nested within a broader function.

meal_sim.f <- function(days, num_hh, ind_hh, iact, group){
  if (group == 'Part-Time') {
    group = "Part-Time"
  } else if (group == "Full-Time") {
    group = "Full-Time"
  } else {
    stop("Argument 'group' must be either Full-Time or Part-Time. Else modify function.")
  }
  if(class(days) == "character"){stop("Error: Days must be numeric.")}
  if(class(num_hh) == "character"){stop("Error: Number of households must be numeric.")}
  if(class(ind_hh)== "character"){stop("Error: Number of individuals per household must be numeric.")}
  if(class(iact) != "data.frame"){stop("Error: 'iact' must be a data frame of Sleep and Eating and Drinking activities.")}
  
  # Requires iact specification from sleep function
  num_hh <- num_hh
  num_ind <- ind_hh
  samplesize <- num_hh * num_ind
  
  ## Set Group Parameters:
  
  if(group == "Part-Time"){
    # Start times
    betasm1 <- c(rate = params$Params[params$Pop=="Part-Time"&params$Var=="Meal 1"&params$Labs=="rate"],
               pi = params$Params[params$Pop=="Part-Time"&params$Var=="Meal 1"&params$Labs=="pi"])
    betasm2 <- c(shape = params$Params[params$Pop=="Part-Time"&params$Var=="Meal 2"&params$Labs=="shape"],
                 rate = params$Params[params$Pop=="Part-Time"&params$Var=="Meal 2"&params$Labs=="rate"])
    
    # Duration
    # Could sample meal durations from observed empirical distributions. Here sampled from observed means.
    mm1_end=rnorm(n=1,40.7,sd=10)
    mm2_end=rnorm(n=1,40.7,sd=10)
    mm3_end=rnorm(n=1,40.7,sd=10)
  }
  
  if(group == "Full-Time"){
    # Start times
    betasm1 <- c(rate = params$Params[params$Pop=="Full-Time"&params$Var=="Meal 1"&params$Labs=="rate"],
                                      pi = params$Params[params$Pop=="Full-Time"&params$Var=="Meal 1"&params$Labs=="pi"])
    betasm2 <- c(shape = params$Params[params$Pop=="Full-Time"&params$Var=="Meal 2"&params$Labs=="shape"],
                 rate = params$Params[params$Pop=="Full-Time"&params$Var=="Meal 2"&params$Labs=="rate"])
    # Durations
    mm1_end=rnorm(n=1,34.4,sd=10)
    mm2_end=rnorm(n=1,34.4,sd=10)
    mm3_end=rnorm(n=1,34.4,sd=10)
  }
  
  # Shared Parameters for Meal 3 transtion
  betasm3 <- c(mean = params$Params[params$Pop=="Full-Time"&params$Var=="Meal 3"&params$Labs=="mean"],
               sd = params$Params[params$Pop=="Full-Time"&params$Var=="Meal 3"&params$Labs=="sd"])
  
  ## Custom Hazard Functions ##
  #Meal 1 - Zero-inflated Exponential
  ZIE_haz.f <- function(x, rate,pi){
    res <- dZIE(x,rate=rate,pi=pi)/pZIE(x,rate=rate,pi=pi,lower.tail = F)
  }
  haz_meal1.f <- function(t,x,betas){
    res <- ZIE_haz.f(x = t,rate=betasm1[1],pi=betasm1[2])
  }
  
  #Meal 2 Gamma
  G_haz.f <- function(x, shape,rate){
    res <- dgamma(x,shape=shape,rate=rate)/pgamma(x,shape=shape,rate=rate,lower.tail = F)
  }
  haz_meal2.f <- function(t,x,betas){
    res <- G_haz.f(x = t,shape=betasm2[1],rate=betasm2[2])
  }
  
  #Meal 3 Gaussian
  N_haz.f <- function(x,mean,sd){
    res <- dnorm(x,mean=mean,sd=sd)/pnorm(x,mean=mean,sd=sd, lower.tail = F)
  }
  haz_meal3.f <- function(t,x,betas){
    res <- N_haz.f(x=t,mean=betasm3[1],sd=betasm3[2])
  }
  
  ## Simulation ## 
  # Simulate over consecutive days:
  temp = data.frame()
  # First meal: #1:(days-1)
  for (i in 1:days){
    # Sleep Timing To Check
    prev_sleep <- iact %>%
      dplyr::filter(event == i-1)%>%
      dplyr::pull(tstop)
    next_sleep <- iact %>% 
      filter(event==i) %>% 
      dplyr::pull(tstart)
    mm1 <- simsurv(betas = betasm1, 
                   x = data.frame(id=1:samplesize), 
                   hazard = haz_meal1.f, 
                   maxt = NULL,
                   interval = c(0,1440*days))
    mm1_start <- ifelse(mm1$eventtime<=1,0,mm1$eventtime)+prev_sleep
    tstopm1 <- mm1_start+mm1_end
    
    # Next Meal
    mm2 <- simsurv(betas = betasm2, 
                     x = data.frame(id=1:samplesize), 
                     hazard = haz_meal2.f, 
                     maxt = NULL,
                     interval = c(0,1440*days))
    mm2_start <- mm2s <- ifelse(mm2$eventtime<=1,0,mm2$eventtime) + tstopm1
    # Fixed Cross-pooling: Force logic to check individual boundaries
    violators2 <- which(mm2_start >= next_sleep)
    if(length(violators2) > 0){
      resampm1 <- simsurv(betas = betasm2,
                              x = data.frame(id=1:100),
                              hazard = haz_meal2.f,
                              maxt = NULL,
                              interval = c(0, 1440*days))
      rs <- as.vector(outer(resampm1$eventtime,tstopm1,"+"))
      for(j in violators2){
        a_idx <- j
        prs <- rs %>%
          subset(rs > tstopm1[a_idx] & rs < next_sleep[a_idx])
        if(length(prs)==1){
          mm2_start[a_idx]<-prs
        } else if(length(prs)>1){
          mm2_start[a_idx] <- sample(prs,size = 1,replace = F)
        } else if(length(prs)==0){
          mm2_start[a_idx] <- next_sleep[a_idx]-mm2_end
        }
        else {stop(paste("Erroring out at PRS functions",j))}
      }
    }
    tstopm2 <- pmin(mm2_start+mm2_end,next_sleep)
    
    # Third Meal
    # No third meal if second meal occurs within one hour of sleep initiation
    meal3ids <- which((next_sleep - tstopm2) > 60)
    # Simulate -- Truncate the meal distribution support by actual available time:
    # Take a block draw, sort, and sample within allowable range.
    mm3 <- simsurv(betas = betasm3,
                   x = data.frame(id=meal3ids),
                   hazard = haz_meal3.f,
                   maxt = NULL,
                   interval = c(0,1440*days))
    mm3_start <- mm3s <- ifelse(mm3$eventtime<=1,0,mm3$eventtime)+tstopm2[meal3ids] 
    #Different censoring rates could be incorporated here
    
    violators3 <- which(mm3_start >= next_sleep[meal3ids])
    if(length(violators3) > 0){
      resampm3 <- simsurv(betas = betasm3,
                          x = data.frame(id=1:100),
                          hazard = haz_meal3.f,
                          maxt = NULL,
                          interval = c(0, 1440*days))
      rs3 <- as.vector(outer(resampm3$eventtime,tstopm2,"+"))
      for(k in violators3){
        actual_idx <- meal3ids[k]
          prs3 <- rs3 %>%
            subset(rs3 > tstopm2[actual_idx] & rs3 < next_sleep[actual_idx])
          if(length(prs3)==1){
            mm3_start[k]<-prs3
          } else if(length(prs3)>1){
            mm3_start[k] <- sample(prs3,size = 1,replace = F)
          } else if(length(prs3)==0){
            mm3_start[k] <- mm3_start-mm3_end
          }
          else {
            stop("Erroring out at PRS3 functions")
          }
      }
    }
      tstopm3 <- pmin(mm3_start + mm3_end, next_sleep[meal3ids]) 
    
    # Bind it all together:
    m_ <- iact %>% 
      dplyr::filter(event==i) %>%
      dplyr::select(id, HH, sex, event)
    mmeal1 <- m_ %>%
      mutate(tstart = mm1_start,
             tstop = tstopm1,
             activity = rep("Eating and Drinking 1",samplesize))
    mmeal2 <- m_ %>%
      mutate(tstart = mm2_start,
             tstop = tstopm2,
             activity = rep("Eating and Drinking 2",samplesize))
    mmeal3 <- m_ %>%
      dplyr::filter(id %in% meal3ids)%>%
      mutate(tstart = mm3_start,
             tstop = tstopm3,
             activity = rep("Eating and Drinking 3",length(meal3ids)))
    # temp <- rbind(temp,mmeal1,mmeal2)%>%
    #   dplyr::filter(tstart != tstop)
    temp <- rbind(temp,mmeal1,mmeal2,mmeal3)%>%
      dplyr::filter(tstart != tstop)
    print(paste("Day ",i,"meals complete; ",days-i," to go."))
  }
  return(temp)
}

# ----------------------
# # Activity Gaps # # 
# ----------------------

# Competing risks to fill in the gaps:
# Not going to use transition probabilities, instead use "period probabilities" 
# that is, activities that occur between events:
# period 1: sleep and meal 1
# period 2: meal 1 and last meal
# period 3: meal 3 and sleep
# All period probabilities based on observed probabilities from empirical time use data.

# --------------------------
# Competing Risks Function #
# --------------------------

# Purpose: generate gap activities (start times, durations) for individuals
# using a competing risk framework. Not currently clustered; could introduce clustering. 
# Ultimately nested within a broader function.

## 

crisk.f <- function(group, iact2, num_hh, ind_hh, days, actlist){
  
  if(!is.data.frame(iact2)){
    stop("Error: 'iact2' must be a data frame of Sleep and Eating and Drinking activities.")
  }
  if(any(stringr::str_detect(actlist, "Sleeping"))){
    stop("Error: list of activities cannot contain sleep activities.")
  }
  if(any(stringr::str_detect(actlist, "Eating"))){
    stop("Error: list of activities cannot contain meals/eating activities.")
  }
  
  num_hh <- num_hh
  num_ind <- ind_hh
  samplesize <- num_hh * num_ind  
  
  ## Key Demographics Mapping
  key <- data.frame(id = 1:samplesize,
                    HH = rep(1:num_hh, each = num_ind),
                    sex = rep(c("M","F"), length.out = samplesize))
  
  ## --- Identify Gaps ----
  iact_2 <- iact2 %>% 
    arrange(id, event, tstart) %>%
    filter(tstop != 0) %>%
    group_by(id) %>%
    mutate(
      Gap_No = paste0("X", row_number()),
      Gap = tstart - lag(tstop),
      Activity = lag(activity),      # Coming from activity...
      Gap_Start = lag(tstop)         # Stop time of previous activity
    ) %>%
    ungroup() %>%
    filter(!is.na(Gap)) %>%          # Drop the first row per person since it has no prior lag
    dplyr::select(id, Gap_No, Gap, Activity, tstart = Gap_Start)
  
  gaps <- iact_2
  
  ## Simulate survival times and activity draws for each period:
  p1 <- p2 <- p3 <- data.frame(time = numeric(), 
                               activity = character(), 
                               id = integer(), 
                               tstart = numeric(), 
                               event = integer())
  
  if(nrow(gaps) > 0) {
    for (i in 1:nrow(gaps)){
      # Extract safe digit ID even if gaps exceed 10+ rows
      current_event <- as.numeric(gsub("\\D", "", gaps$Gap_No[i])) - 1
      
      # --- PERIOD 1: Sleeping ---
      if(gaps$Activity[i] == "Sleeping"){
        # Get Hazard Rates for Groups:
        if (group == "Part-Time"){
          pmat <- dat$PTP1
        } else if (group == "Full-Time") {
          pmat <- dat$FTP1
        } else {
          stop("Invalid group designation.")
        }
        
        # Initialize
        m <- 0
        t <- c()
        locs <- c()
        events <- data.frame()
        while (m < gaps$Gap[i]){
          events <- rexp(dim(pmat)[2],pmat)*gaps$Gap[i]
          # for(k in 1:dim(pmat)[2]){
          #   e <- simsurv(dist="exponential",
          #                        lambda = pmat[k],
          #                        x = data.frame(id=1),
          #                interval = c(1, gaps$Gap[i]))
          #   events <- rbind(events,data.frame(e))%>%
          #     mutate(eventtime = ifelse(eventtime < 1, 1, eventtime))
          # }
          loc <- which(events == min(events))[1]
          event.time <- min(events,na.rm=T)
          m = m + event.time
          t = append(t, event.time)
          locs <- append(locs,loc)
        }
        if(!is.null(t) && !is.na(t[1]) && t[1] <= gaps$Gap[i]){
          # Isolate pieces safely matching length requirements
          sub_times <- c(t[1:(length(t)-1)], gaps$Gap[i] - sum(t[1:(length(t)-1)], na.rm = TRUE))
          
          temp1 = data.frame(time = sub_times,
                             activity = labels(pmat)[[2]][locs],
                             id = rep(gaps$id[i], length(sub_times)),
                             event = current_event)
          
          # Reconstruct step intervals incrementally
          temp1$tstart <- gaps$tstart[i] + cumsum(c(0, sub_times[-length(sub_times)]))
          p1 = rbind(p1, temp1)
        }
        if(!is.null(t) && t[1] > gaps$Gap[i] || !is.null(t) && is.na(t[1])){
          temp1 = data.frame(time = gaps$Gap[i],
                             activity = labels(pmat)[[2]][locs][1],
                             id = gaps$id[i],
                             event = current_event,
                             tstart = gaps$tstart[i])
          p1 = rbind(p1, temp1)
        }
        if(is.null(t)){
          temp1 = data.frame(time = gaps$Gap[i],
                             activity = gaps$Activity[i],
                             id = gaps$id[i],
                             event = current_event,
                             tstart = gaps$tstart[i])
          p1 = rbind(p1, temp1)
        }
      }
      
      # --- PERIOD 3: Eating and Drinking 3 ---
      if(gaps$Activity[i] == "Eating and Drinking 3"){
        if (group == "Part-Time"){
          pmat <- dat$PTP3[dat$PTP3 > 0]
          names(pmat) <- c(labels(dat$PTP3)[[2]][1:length(pmat)])
        } else if (group == "Full-Time") {
          pmat <- dat$FTP3[dat$FTP3 > 0]
          names(pmat) <- c(labels(dat$FTP3)[[2]][1:length(pmat)])
        } else {
          stop("Invalid group designation.")
        }
        m3 <- 0
        t3 <- c()
        locs3 <- c()
        while (m3 < gaps$Gap[i]){
          events3 <- rexp(length(pmat),pmat)*gaps$Gap[i]
          loc3 <- which(events3 == min(events3))[1]
          event.time3 <- min(events3,na.rm=T)
          m3 = m3 + event.time3
          t3 = append(t3, event.time3)
          locs3 <- append(locs3,loc3)
        }
        if(!is.null(t3) && !is.na(t3[2]) && t3[2] <= gaps$Gap[i] ){
          sub_times3 <- c(t3[1:(length(t3)-1)], gaps$Gap[i] - sum(t3[1:(length(t3)-1)], na.rm=T))
          
          temp3 = data.frame(time = sub_times3,
                             activity = names(pmat)[locs3],
                             id = rep(gaps$id[i], length(sub_times3)),
                             event = current_event)
          temp3$tstart <- gaps$tstart[i] + cumsum(c(0, sub_times3[-length(sub_times3)]))
          p3 = rbind(p3, temp3)
        }
        if(!is.null(t3) && t3[2] > gaps$Gap[i] || !is.null(t3) && is.na(t3[2])){
          temp3 = data.frame(time = gaps$Gap[i],
                             activity = names(pmat)[locs3][1],
                             id = gaps$id[i],
                             event = current_event,
                             tstart = gaps$tstart[i])
          p3 = rbind(p3, temp3)
        }
        if(is.null(t)){
          temp3 = data.frame(time = gaps$Gap[i],
                             activity = gaps$Activity[i],
                             id = gaps$id[i],
                             event = current_event,
                             tstart = gaps$tstart[i])
          p3 = rbind(p3, temp3)
        }
      }
      
      # --- PERIOD 2: Eating and Drinking 1 or 2 ---
      if(gaps$Activity[i] == "Eating and Drinking 1" || gaps$Activity[i] == "Eating and Drinking 2") {
        if (group == "Part-Time"){
          pmat <- dat$PTP2
        } else if (group == "Full-Time") {
          pmat <- dat$FTP2
        } else {
          stop("Invalid group designation.")
        }
        m2 <- 0
        t2 <- c()
        locs2 <- c()
        while (m2 < gaps$Gap[i]){
          events2 <- rexp(dim(pmat)[2],pmat)*gaps$Gap[i]
          # for(o in 1:dim(pmat)[2]){
          #   e2 <- simsurv(dist="exponential",
          #                lambda = pmat[o],
          #                x = data.frame(id=1),
          #                interval=c(1,gaps$Gap[i]))
          #   events2 <- rbind(events2,data.frame(e2))%>%
          #     mutate(eventtime = ifelse(eventtime < 1, 1, eventtime))
          # }
          loc2 <- which(events2 == min(events2))[1]
          event.time2 <- min(events2,na.rm=T)
          locs2 <- append(locs2,loc2)
          m2 = m2+event.time2+1
          t2 = append(t2, event.time2)
        }
        if(!is.null(t2) && !is.na(t2[2]) && t2[2] <= gaps$Gap[i] ){
          sub_times2 <- c(t2[1:(length(t2)-1)], gaps$Gap[i] - sum(t2[1:(length(t2)-1)], na.rm = TRUE))
          
          temp2 = data.frame(time = sub_times2,
                             activity = labels(pmat)[[2]][locs2],
                             id = rep(gaps$id[i], length(sub_times2)),
                             event = current_event)
          temp2$tstart <- gaps$tstart[i] + cumsum(c(0, sub_times2[-length(sub_times2)]))
          p2 = rbind(p2, temp2)
        }
        if(!is.null(t2) && t2[2] > gaps$Gap[i] || !is.null(t2) && is.na(t2[2])){
          temp2 = data.frame(time = gaps$Gap[i],
                             activity = labels(pmat)[[2]][locs2][1],
                             id = gaps$id[i],
                             event = current_event,
                             tstart = gaps$tstart[i])
          p2 = rbind(p2, temp2)
        }
        if(is.null(t2)){
          temp2 = data.frame(time = gaps$Gap[i],
                             activity = gaps$Activity[i],
                             id = gaps$id[i],
                             event = current_event,
                             tstart = gaps$tstart[i])
          p2 = rbind(p2, temp2)
        }
      }
      #if(floor(i/nrow(gaps)) %% .25 == 0){print(paste(round(i/nrow(gaps),2)*100,"% gaps filled.",sep=""))} 
    }
  }
  
  # Bind arrays, explicitly structure tstop column data, drop scratch variables
  res <- rbind(p1, p2, p3) %>%
    mutate(tstop = tstart + time) %>%
    dplyr::select(id, tstart, tstop, activity, event)%>%
    dplyr::filter(tstart != tstop)
  
  # Append baseline household structure demographic details
  res <- res %>% left_join(key, by = "id")
  row.names(res) <- NULL
  
  return(res)
}

# ------------------------------
# # Full Simulation Function # #
# ------------------------------

big_sim.f <- function(days, num_hh, ind_hh, group,actlist,params){
  message("HASP has begun. Your time use simulations should be ready shortly.")
  iact <- sleep_sim.f(days = days, num_hh=num_hh, ind_hh=ind_hh, group=group,params=params)
  iact2 <- meal_sim.f(days = days, num_hh=num_hh, ind_hh=ind_hh,iact = iact, group=group)
  iact2 <- rbind(iact,iact2)
  iact3 <- crisk.f(group = group, iact2, num_hh = num_hh, ind_hh = ind_hh, days = days,actlist)
  iact_t <- rbind(iact2,iact3) %>% 
    filter(tstart != 0 & tstop != 0) %>%
    filter(tstart < tstop) # Filter out last activity line (not a full day)
  return(iact_t)
}

# ------------------------------
# # Compute Details Function # #
# ------------------------------

compute_details.f <- function(func, ...) { 
  message("Analyzing target function profile (Single-Pass)... Please wait.") 
  
  # 1. Structural Details (Static Code Analysis) 
  func_body <- body(func) 
  body_lines <- deparse(func_body) 
  line_count <- length(body_lines) 
  loop_count <- sum(stringr::str_count(body_lines, "\\b(for|while)\\b")) 
  
  # 2. Performance & Memory Profiling
  # Clear existing garbage and reset the high-water mark counters
  gc(verbose = FALSE, full = TRUE)
  gc(reset = TRUE) 
  
  # Execute the target function exactly ONCE and record time 
  time_taken <- system.time({ 
    result <- func(...) 
  }) 
  
  # Read the peak memory stats immediately after running
  mem_post <- gc()
  
  # Row 1 is Ncells, Row 2 is Vcells. Column 4 is "max used (Mb)"
  peak_ncells_mb <- mem_post[1, 4]
  peak_vcells_mb <- mem_post[2, 4]
  total_peak_mb  <- round(peak_ncells_mb + peak_vcells_mb, 2)
  
  # 3. Compile Meta Metrics 
  meta <- list( 
    structural = list( 
      line_count = line_count, 
      loops_detected = loop_count 
    ), 
    performance = list( 
      execution_time_ms = as.numeric(time_taken["elapsed"]) * 1000, 
      approx_peak_mem_mb = paste0(total_peak_mb, " MB") 
    ) 
  ) 
  
  # Return metrics along with the original output safely 
  return(list(metrics = meta, output = result)) 
}

# --------------------
# # Data Cleaning # #
# -------------------

# Purpose: Separate simulated continuous time activities into discrete days (day sect),
# mostly for plotting an analysis.

daysect.f <- function(df){
  ## Isolate problem rows
  df_issue <- df %>%
    filter(floor(tstart/1440) != floor(tstop/1440))%>%
    mutate(tstart1 = tstart,
           tstop1 = floor(tstop/1440)*1440,
           tstart2 = floor(tstop/1440)*1440,
           tstop2 = tstop) %>%
    dplyr::select(id,HH,sex,activity,event,tstart1,tstart2,tstop1,tstop2)
  ## Get previous day start/stop
  prior <- df_issue %>%
    dplyr::select(-starts_with("tstop"))%>% 
    pivot_longer(cols = -c(id,HH,sex,activity,event),
                 names_to = "type",
                 values_to = "tstart")%>%
    dplyr::select(-type)%>%
    arrange(id,HH,sex,activity,event)
  ## Get next day start/stop
  post <- df_issue %>%
    dplyr::select(-starts_with("tstart"))%>%
    pivot_longer(cols = -c(id,HH,sex,activity,event),
                 names_to = "type",
                 values_to = "tstop")%>%
    dplyr::select(-type)%>%
    arrange(id,HH,sex,activity,event)
  ## Merge total
  total <- cbind(prior,tstop = post$tstop)
  ## Re collate
  df_clean <- df %>%
    filter(floor(tstart/1440) == floor(tstop/1440))%>%
    rbind(., total)%>%
    arrange(HH,id,sex,event,tstart)
}



## OPEN ISSUES
# 1)  Weighting of ATE in model where an individual ATE is calculated for each 
#     treated unit (right now, using average pre-treatment weight)
# 2)  Covariates


## Preliminaries

## Set up libraries
library(tidyverse)
library(boot)
library(data.table)
library(synthdid)
library(modelsummary)
library(fixest)
library(huxtable)

## Set up directories - change to local directory
directory <- "/Users/johniselin/Library/CloudStorage/Dropbox/sdid_weights"
results <- file.path(directory,"results")
data <- file.path(directory,"data")

setwd(directory)

## Set Seed 
set.seed(5)

# Define Monte Carlo simulation parameters
N_states <- 100  # Number of states
T_periods <- 20   # Number of time periods
N_simulations <- 100  # Number of simulations
treatment_period <- 10  # Time period after which treatment is allocated
treatment_state <- 20  # Max state number for treatment (i.e. states 1 - X are treated)
treatment_size <- 2 # True treatment effect 
bootstrap_size <- 50

### DEFINE FUNCTIONS 

## DATA GENERATION 
f_data_gen <- function(N_states, 
                       T_periods, 
                       treatment_period, 
                       treatment_type = 1, 
                       treatment_size = 1, 
                       state_gen_struc = 1, 
                       state_fe_m = 5, 
                       state_fe_sd = 1,
                       time_fe_m = 3, 
                       time_fe_sd = 1,
                       pop_min = 10,
                       pop_max = 20, 
                       rate = 0.02) {
  
  # Generate synthetic panel data
  tmp_data <- expand.grid(state = 1:N_states,       # Number of States
                      period = 1:T_periods) %>% # Number of Periods 
    mutate(
      # Assign treatment status to states
      treated = ifelse(state <= treatment_state,1,0),
      # Assign treatment status to state X period  
      did = ifelse(state <= treatment_state & period >= treatment_period, 1, 0), 
      # Generate random error term
      e = rnorm(n(), mean = 0, sd = 1)
    ) %>%
    ## Group by state and create state FE
    group_by(state) %>% 
    mutate(state_fe = rnorm(1, mean = state_fe_m, sd = state_fe_m)) %>%
    ungroup() %>%
    ## Group by period and create time FE
    group_by(period) %>% 
    mutate(period_fe = rnorm(1, mean = time_fe_m, sd = time_fe_sd)) %>%
    ungroup() %>%
    ## Group by state and create population 
    arrange(state, period) %>%
    group_by(state) %>%
    mutate(tmp = runif(1,0,1), 
           pop_initial = runif(1,pop_min, pop_max)) %>%
    mutate(r = ifelse(tmp > 0.5, rate, rate * (-1))) %>%
    mutate(pop = pop_initial * (1 + r)^period) %>%
    ungroup() %>%
    select(-tmp, -pop_initial, -r)  
  
    
    # Homogeneous treatment - same effect size across state and period
    if (treatment_type == 1) {
      
      tmp_data <- tmp_data %>% mutate(treatment_effect = did * treatment_size) 
      
    }
  
    # Heterogeneous treatment - Effect varies with state FE 
    else if (treatment_type == 2) {
    
      tmp_data <- tmp_data %>%
        mutate(treatment_effect = did * treatment_size * (1 + (state_fe-state_fe_m)/state_fe_m))
    }    
  
    # Heterogeneous treatment - Effect varies with state weight 
    else if (treatment_type == 3) {
    
      tmp_data <- tmp_data %>%
        group_by(period) %>% mutate(pop_mean = mean(pop)) %>% ungroup() %>%
        mutate(
          treatment_effect = did * treatment_size * (1 + (pop - pop_mean)/pop_mean)) %>%
        select(-pop_mean) 
    
    }        

    # Additive FE 
    if (state_gen_struc == 1) {
  
      tmp_data <- tmp_data %>%
        mutate( y = state_fe + period_fe + e + treatment_effect)
    
    }  
      
    # Interactive FE
    else if (state_gen_struc == 2) {
    
      tmp_data <- tmp_data %>%
        mutate( y = state_fe + period_fe + state_fe * period_fe + e + treatment_effect)

    }  
  
  # Combine
  return(tmp_data)
  
}

## DEFINE COVARIATEs SDID SETUP FUNCTION 
setup_function <- function(data, 
                           unit, 
                           time, 
                           treatment, 
                           covar_list){
  
  table <- data 
  
  colnames(table)[colnames(table) == unit] <- 'unit'
  colnames(table)[colnames(table) == time] <- 'time'
  colnames(table)[colnames(table) == outcome] <- 'outcome'
  colnames(table)[colnames(table) == treatment] <- 'treatment'
  
  table <- table %>% group_by(unit) %>% mutate(ever_treat = max(treatment))
  table <- table %>% arrange(ever_treat, unit, time)
  
  Y <- array(table[['outcome']], dim=c(N_states, T_periods, 1))
  N0 <- table %>% group_by(unit) %>% summarize(ever_treat = max(ever_treat)) %>% ungroup() %>% summarize(count = sum(ever_treat)) %>% as.numeric()
  T0 <- treatment_period - 1
  
  table <- table %>% arrange(treatment, time, unit)
  
  list_covs <- c()
  for(cv in covar_list){
    list_covs <- c(list_covs, as.vector(table[[cv]]))
  }
  
  ## Create new version of table sorted by treatment, time, and unit 
  table <- data %>% arrange(treatment, time, unit)
  
  ## Create empty list 
  list_covs <- c()
  
  ## Loop over each covariate 
  for(cv in covar_list){
    
    ## Add covariate to list 
    list_covs <- c(list_covs, as.vector(table[[cv]]))
    
  }
  
  ## Create array of correct format 
  x <- array(list_covs, dim=c(N_states, T_periods, length(covar_list)))
  
  return(x)
}


## DEFINE BOOTSTRAP FUNCTION FOR WEIGHTED AVERAGE VERSION 
boot_model1 <- function(data, indices ){
  
  ## Set up data (in wide format so that rows = states)
  data <- data[indices,]
  
  ## Reshape long (state X period)
  reshape_data <- data %>%
    pivot_longer(cols = !state, 
                 names_to = c(".value", "period"), 
                 names_pattern = "(.*)_(.*)") %>%
    data.frame(reshape_data) %>%
    mutate(period = as.numeric(period))
  
  ## Check that there are treated and control units 
  temp3 <- reshape_data %>% group_by(treated) %>% summarize(count = n())
  t <- as.numeric(temp3[1,2])
  c <- as.numeric(temp3[2,2])
  
  if(t == 0 | c == 0 ) {
    print("Skipped")
  }
  else {
    
    ## Take the weighted average of the treated units
    avg_data <- reshape_data %>%
      mutate(weight = !! sym(weight_var)) %>%
      mutate(collapse_var = ifelse(treated == 1, 1, state)) %>%
      group_by(collapse_var, period) %>%
      summarise(y = weighted.mean(y,weight), 
                did = mean(did), 
                treated = mean(treated),
                .groups = "drop_last") %>%
      ungroup() %>%
      rename(state = collapse_var)
    
    ## Basic set-up 
    setup <- avg_data %>%
      dplyr::select(state, period, y, did) %>%
      as.data.frame() %>%
      synthdid::panel.matrices(
        unit = 1,
        time = 2,
        outcome = 3,
        treatment = 4
      )
    
    ## RUN SDID 
    est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
    sum_est <- summary(est)
    ate <- sum_est$estimate
    return(ate)
    
  }
  
}


## DEFINE BOOTSTRAP FUNCTION FOR WEIGHTED AVERAGE OF INDIVIDUAL SDID VERSION
boot_model2 <- function(data, indices ){
  
  ## Set up data (in wide format so that rows = states)
  data <- data[indices,]
  
  ## Reshape long (state X period)
  reshape_data <- data %>%
    pivot_longer(cols = !state, 
                 names_to = c(".value", "period"), 
                 names_pattern = "(.*)_(.*)") %>%
    data.frame(reshape_data) %>%
    mutate(period = as.numeric(period))
  
  ## Check that there are treated and control units 
  temp3 <- reshape_data %>% group_by(treated) %>% summarize(count = n())
  t <- as.numeric(temp3[1,2])
  c <- as.numeric(temp3[2,2])
  
  if(t == 0 | c == 0 ) {
    print("Skipped")
  }
  else {
    
    ## Construct a dataset of treated units and their weights
    treated <- data %>% 
      mutate(weight = !!sym(weight_var)) %>%
      filter(treated == 1 & did == 0) %>% 
      select(state, period, weight) %>% 
      group_by(state) %>%
      summarise(weight = mean(weight)) %>%
      ungroup() %>%
      mutate(ate = NA)
    
    ## Construct a list of treated states
    treated_list <- treated %>% select(state) %>% unique() %>% as.list()
    treated_list <- treated_list[[1]]
    
    ## Loop over each treated state 
    for (s in treated_list) {
      
      ## Basic set-up, keeping only one treated state + control states
      setup <- data %>%
        filter(treated == 0 | state == s) %>%
        dplyr::select(state, period, y, did) %>%
        as.data.frame() %>%
        synthdid::panel.matrices(
          unit = 1,
          time = 2,
          outcome = 3,
          treatment = 4
        )
      
      ## Run SDID 
      tmp_est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
      tmp_sum_est <- summary(tmp_est)
      tmp_ate <- tmp_sum_est$estimate
      
      ## Add ATE to dataset 
      treated <- treated %>% mutate(ate = if_else(state == s, tmp_ate, ate))
      
    }
    
    ## Get weighted average of ATE 
    ate <- treated %>%
      summarise(ate = weighted.mean(ate,weight)) %>%
      ungroup() %>%
      as.numeric()

    return(ate)
    
  }
  
}


## DEFINE SDID FUNCTION 
f_sdid <- function(data, 
                   weighted = 0, 
                   weight_var = FALSE, 
                   covar_list = FALSE, # MAKE EMPTY LIST
                   calc_se = 0, # 0 = no SEs, 1 = calculate SEs
                   n_boot = 100
                   ) {

  ## OPTION 1: No weights   
  if (weighted == 0) { 
    
    ## Display weights 
    print("Run without weights")
    
    ## Basic set-up 
   setup <- data %>%
     dplyr::select(state, period, y, did) %>%
     as.data.frame() %>%
     synthdid::panel.matrices(
       unit = 1,
       time = 2,
       outcome = 3,
       treatment = 4
      )
   
    ## COVARIATES 
    if (covar_list != FALSE) {
      
      ## Create array of covariates 
      x <- setup_function(data, state, period, treated, covar_list = covar_list)
      
      ## Run SDID 
      est <- synthdid_estimate(setup$Y, setup$N0, setup$T0, X = x)     
    }
      
    ## NO COVARIATES 
    else {
    
      ## Run SDID 
       est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)  
       
    }
   
    ## Run SDID 
    sum_est <- summary(est)
       
    ## Store ATE 
    ate <- sum_est$estimate
    
    ## Store SE 
    if (calc_se == 1) {
      se <- sum_est$se
    }
    else {
      se = NA
    } 
    
    ## Store list for export
    val <- c(ate, se)
  }
  
  ## OPTION 2: Run using weights, using a weighted average treated unit 
  if (weighted == 1) { 
    
    ## Check that weighting variable was specified
    if (weight_var == FALSE) {
      print("No weighting variable specified")
      break 
    }
    else {
      print(paste("Weighting variable ", weight_var))
    }
    
    ## Take the weighted average of the treated units 
    avg_data <- data %>%
      mutate(collapse_var = ifelse(treated == 1, 1, state)) %>%
      mutate(weight = !! sym(weight_var)) %>%
      group_by(collapse_var, period) %>%
      summarise(y = weighted.mean(y, weight), 
                did = mean(did), 
                treated = mean(treated), 
                .groups = "keep") %>%
      ungroup() %>%
      rename(state = collapse_var)
    
    ## Basic set-up 
    setup <- avg_data %>%
      dplyr::select(state, period, y, did) %>%
      as.data.frame() %>%
      synthdid::panel.matrices(
        unit = 1,
        time = 2,
        outcome = 3,
        treatment = 4
      )
    
    ## COVARIATES 
    if (covar_list != FALSE) {
      
      ## Create array of covariates 
      x <- setup_function(data, state, period, treated, covar_list = covar_list)
      
      ## Run SDID 
      est <- synthdid_estimate(setup$Y, setup$N0, setup$T0, X = x)     
    }
    
    ## NO COVARIATES 
    else {
      
      ## Run SDID 
      est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)  
      
    }    
    
    ## RUN SDID 
    est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
    
    
    sum_est <- summary(est)
    
    ## Store ATE 
    ate <- sum_est$estimate
    
    ## Store SE 
    if (calc_se == 1) {
      
      ## Reshape data so that the boot function index option works 
      boot_data <- data %>%
        pivot_wider(id_cols = state, names_from = period, values_from = c(y, pop, did, treated))
      
      ## Run boostrap command
      boot_est <- boot(boot_data, boot_model1, R = n_boot, stype = "i")
     
      ## Store SE
      se <- sd(boot_est$t)
      
    }
    else {
      se = NA
    } 
    
    ## Store list for export
    val <- c(ate, se)
  }
    
  ## OPTION 3: Run using weights, constructing a separate SDID model for each treated unit 
  if (weighted == 2) { 
    
    ## Check that weighting variable was specified
    if (weight_var == FALSE) {
      print("No weighting variable specified")
      break 
    }
    else {
      print(paste("Weighting variable ", weight_var))
    }
    
    ## Construct a dataset of treated units and their weights
    treated <- data %>% 
      mutate(weight = !! sym(weight_var)) %>%
      filter(treated == 1 & did == 0) %>% 
      select(state, period, weight) %>% 
      group_by(state) %>%
      summarise(weight = mean(weight)) %>%
      ungroup() %>%
      mutate(ate = NA)
    
    ## Construct a list of treated states
    treated_list <- treated %>% select(state) %>% unique() %>% as.list()
    treated_list <- treated_list[[1]]
    
    ## Loop over each treated state 
    for (s in treated_list) {

      ## Basic set-up, keeping only one treated state + control states
      setup <- data %>%
        filter(treated == 0 | state == s) %>%
        dplyr::select(state, period, y, did) %>%
        as.data.frame() %>%
        synthdid::panel.matrices(
          unit = 1,
          time = 2,
          outcome = 3,
          treatment = 4
        )

      ## Run SDID 
      tmp_est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
      tmp_sum_est <- summary(tmp_est)
      tmp_ate <- tmp_sum_est$estimate

      ## Add ATE to dataset 
      treated <- treated %>% mutate(ate = if_else(state == s, tmp_ate, ate))

          }
    
    ## Get weighted average of ATE 
    ate <- treated %>%
      summarise(ate = weighted.mean(ate,weight)) %>%
      as.numeric()

    
    ## Store SE 
    if (calc_se == 1) {
      
      ## Reshape data so that the boot function index option works 
      boot_data <- data %>%
        mutate(weight = !! sym(weight_var)) %>%
        pivot_wider(id_cols = state, names_from = period, values_from = c(y, weight, did, treated))
      
      ## Run boostrap command
      boot_est <- boot(boot_data, boot_model2, R = n_boot, stype = "i")
      
      ## Store SE
      se <- sd(boot_est$t)
      
    }
    else {
      se = NA
    } 
    
    ## Store both values 
    val <- c(ate, se)
   
  } 
    
    ## RETURN 
    return(val)
    
  }


### RUN SIMULATIONS

## Simulation 1: Homogeneous treatment  

# Generate synthetic panel data
simulation <- expand.grid(sim = 1:N_simulations, treat_type = 1:3, data_type = 1:2) %>%
  mutate(ate = NA, ate_nw = NA, ate_wt_1 = NA, ate_wt_2 = NA, did_nw = NA, did_wt = NA)

## LOOP OVER SIMULATIONS
for (i in 1:N_simulations) {
  
  ## LOOP OVER TREATMENT TYPES
  for (j in 1:3) {
    
    ## LOOP OVER DATA GENERATING TYPES 
    for (k in 1:2) {
      
      print(paste("Simulation", i, j, k))

      ## Generate data 
      df <- f_data_gen(N_states = N_states, 
                       T_periods = T_periods, 
                       treatment_period = treatment_period, 
                       treatment_size = treatment_size, 
                       state_gen_struc = k, 
                       treatment_type = j)
      

      ## Get true ATE
      tmp_ate <- df %>% filter(did == 1) %>% summarise(ate = weighted.mean(treatment_effect, pop))
      tmp_ate <- tmp_ate[[1]]
      
      ## Run unweighted SDID 
      tmp <- f_sdid(df, weighted = 0)
      tmp_ate_nw <- tmp[1]
      
      ## Run weighted SDID (version 1)
      tmp <- f_sdid(df, weighted = 1, weight_var = "pop")
      tmp_ate_wt_1 <- tmp[1]
      
      ## Run weighted SDID (version 2)
      tmp <- f_sdid(df, weighted = 2, weight_var = "pop")
      tmp_ate_wt_2 <- tmp[1]
      
      ## Run DiD without weight 
      tmp <- feols(y ~ did | state + period,
                   data = df)
      tmp_did_nw <- tmp$coefficients
      
      ## Run DiD with weight 
      tmp <- feols(y ~ did | state + period,
                            data = df, weights = ~pop)
      tmp_did_wt <- tmp$coefficients

      simulation <- simulation %>% 
        mutate(ate = if_else(sim == i & treat_type == j & data_type == k, tmp_ate, ate), 
               ate_nw = if_else(sim == i & treat_type == j & data_type == k, tmp_ate_nw, ate_nw), 
               ate_wt_1 = if_else(sim == i & treat_type == j & data_type == k, tmp_ate_wt_1, ate_wt_1), 
               ate_wt_2 = if_else(sim == i & treat_type == j & data_type == k, tmp_ate_wt_2, ate_wt_2),
               did_nw = if_else(sim == i & treat_type == j & data_type == k, tmp_did_nw, did_nw), 
               did_wt = if_else(sim == i & treat_type == j & data_type == k, tmp_did_wt, did_wt)
               )
      
      rm(tmp, tmp_did_nw, tmp_did_wt, tmp_ate, tmp_ate_nw, tmp_ate_wt_1, tmp_ate_wt_2)
    }
  }
}


###############

## Export 
write.csv(simulation, file = paste(data, "data.csv", sep = "/"))






## REPLICATION AND EXTENSION OF BORGSCHULTE AND VOGLER (2020)
## EMPRICAL EXAMPLE OF Role of weighted SDID

## Preliminaries

##devtools::install_github("synth-inference/synthdid")

## Set up libraries
library(tidyverse)
library(boot)
library(data.table)
library(synthdid)
library(modelsummary)
library(fixest)

## Set up directories - change to local directory
directory <- "/Users/johniselin/Library/CloudStorage/Dropbox/sdid_weights"
results <- file.path(directory,"results")
data <- file.path(directory,"data")

setwd(directory)

## Set Seed 
set.seed(5)

## ACA Determination 

## ACA Expansion states (First half of 2014)
aca <- c(4, 5, 6, 8, 9, 10, 11,	15, 17, 19, 21, 24, 25, 26, 27, 32, 34, 
         35, 36, 38, 39, 41, 44, 50, 53, 54)	

## DROP AK IN, LA, MT, NH, and PA for later expansions (post first half of 2014)
drop <- c(2, 18, 22, 30, 33, 42)

## IF 2017 - 2019 period included, drop ME, VA 
drop_2019 <- c(23, 51)

## LOAD DATA 

## Mortality data via https://wonder.cdc.gov/ucd-icd10.html

df_mortality <- read.csv(file = file.path(data,"mortality/cdc_wonder_data.csv")) %>%
  ## FLAG SUPPRESSED VALUES 
  mutate(death_supp = ifelse(Deaths == "Suppressed", 1,0),
         death_miss = ifelse(Deaths == "Missing", 1, 0), 
         pop_miss = ifelse(Population == "Missing", 1, 0)) %>%
  mutate(Deaths = ifelse(Deaths == "Suppressed", 0, Deaths),
         Deaths = ifelse(Deaths == "Missing", 0,Deaths), 
         Population = ifelse(Population == "Missing", 0, Population) ) %>%
  mutate(Deaths = as.numeric(Deaths), Population = as.numeric(Population)) %>%
  mutate(crude_rate = (Deaths / Population) * 1000 ) %>%
  select(-"Crude.Rate") %>%
  rename(state_fips = State.Code, county_fips = County.Code, state_name = State)

## Rename all to lower case 
names(df_mortality) <- tolower(names(df_mortality))

## Population data via https://seer.cancer.gov/popdata/download.html
## 1990-2020, 4 Expanded Races by Origin, All US 
df_pop <- read_fwf(file = file.path(data,"covariates/us.1990_2020.19ages.adjusted.txt"),
                   fwf_widths( c(4, 2, 2, 3, 2, 1, 1, 1, 2, 8),
                               c("year", "state_abb", "state_fips", "county_fips", 
                                 "registry", "race", "hispanic", "sex", "age", "pop"))) %>%
  filter(year >= 2007) %>% filter(year <= 2019) %>%
  mutate(age = as.numeric(age), pop = as.numeric(pop)) %>%
  select(-registry, -hispanic) %>%
  group_by(year, state_abb, state_fips, county_fips, race, sex, age) %>%
  summarize(pop = sum(pop)) %>% group_by() %>%
  mutate(county_fips = as.numeric(paste0(state_fips, county_fips)))

## Create relevant population samples 

## Percent of population that is white 
df_pop_race <- df_pop %>%
  select(year, state_abb, state_fips, county_fips, race, pop) %>%
  group_by(year, state_abb, state_fips, county_fips, race) %>%
  summarize(pop = sum(pop)) %>% group_by() %>%
  pivot_wider(names_from = race, values_from = pop, values_fill = 0) %>%
  rename(white = '1', black = '2', aian = '3', asian = '4') %>%
  replace_na(list(white = 0, black = 0, aian = 0, asian = 0)) %>%
  mutate(pct_white = white / (white + black + aian + asian)) %>%
  select(year, state_abb, state_fips, county_fips, pct_white) 

## By-age statistics
df_pop_age <- df_pop %>%
  select(year, state_abb, state_fips, county_fips, age, pop) %>%
  group_by(year, state_abb, state_fips, county_fips, age) %>%
  summarize(pop = sum(pop)) %>% group_by() %>%
  pivot_wider(names_from = age, values_from = pop, values_fill = 0) %>%
  rename(pop_00 = '0', pop_01_04 = '1', pop_05_09 = '2', pop_10_14 = '3', 
         pop_15_19 = '4', pop_20_24 = '5', pop_25_29 = '6', pop_30_34 = '7',
         pop_35_39 = '8', pop_40_44 = '9', pop_45_49 = '10', pop_50_54 = '11',
         pop_55_59 = '12', pop_60_64 = '13', pop_65_69 = '14', pop_70_74 = '15',
         pop_75_79 = '16', pop_80_84 = '17', pop_85 = '18' ) %>%
  mutate(pop_total = rowSums(across(pop_00:pop_85)), 
         pop_20_64 = rowSums(across(pop_20_24:pop_20_24))) %>%
  mutate(pct_55_64 = (pop_55_59 + pop_60_64)/ pop_total, 
         log_20_64 = log(rowSums(across(pop_20_24:pop_20_24))), 
         log_35_44 = log(rowSums(across(pop_35_39:pop_40_44)))) %>%
  select(year, state_abb, state_fips, county_fips, pct_55_64, log_20_64, log_35_44, pop_20_64, pop_total)

## By-age and gender statistics
df_pop_f_age <- df_pop %>%
  filter(sex == 2) %>%
  select(year, state_abb, state_fips, county_fips, age, pop) %>%
  group_by(year, state_abb, state_fips, county_fips, age) %>%
  summarize(pop = sum(pop)) %>% group_by() %>%
  pivot_wider(names_from = age, values_from = pop, values_fill = 0) %>%
  rename(pop_00 = '0', pop_01_04 = '1', pop_05_09 = '2', pop_10_14 = '3', 
         pop_15_19 = '4', pop_20_24 = '5', pop_25_29 = '6', pop_30_34 = '7',
         pop_35_39 = '8', pop_40_44 = '9', pop_45_49 = '10', pop_50_54 = '11',
         pop_55_59 = '12', pop_60_64 = '13', pop_65_69 = '14', pop_70_74 = '15',
         pop_75_79 = '16', pop_80_84 = '17', pop_85 = '18' ) %>%
  mutate(pop_total = rowSums(across(pop_00:pop_85))) %>%
  mutate(log_f_20_64 = log(rowSums(across(pop_20_24:pop_60_64)))) %>%
  select(year, state_abb, state_fips, county_fips, log_f_20_64)

## Merge together population covariates 
df_pop_cov <- df_pop_race %>% 
  full_join(df_pop_age, 
            by = join_by(year, state_abb, state_fips, county_fips)) %>% 
  left_join(df_pop_f_age, 
            by = join_by(year, state_abb, state_fips, county_fips)) %>% 
  mutate(state_fips = as.numeric(state_fips))

## Unemployment data via https://download.bls.gov/pub/time.series/la/
## la.data.64.County
df_unempl <- read_tsv(file = file.path(data,"covariates/la.data.64.County")) %>%
  filter(year >= 2007, year <= 2019) %>%
  select(-footnote_codes) %>%
  mutate(series = as.numeric(substr(series_id, 19, 20)), 
         state_fips = substr(series_id, 6, 7),
         county_fips = substr(series_id, 8, 10)) %>%
  filter(series == 3) %>%
  group_by(state_fips, county_fips, series_id, year) %>%
  summarise(unemp = mean(value)) %>%
  group_by() %>%
  select(year, state_fips, county_fips, unemp) %>%
  mutate(county_fips = as.numeric(paste0(state_fips, county_fips)), 
         state_fips = as.numeric(state_fips))

## Combine data 
df <- df_pop_cov %>% 
  left_join(df_unempl, join_by(year, state_fips, county_fips)) %>%
  left_join(df_mortality, join_by(year, state_fips, county_fips))

## Export version of data 
write.csv(df,file = file.path(data,"analysis_data.csv"))



## Remove intermediate datasets 
rm(df_pop, df_pop_age, df_pop_cov, df_pop_f_age, df_pop_race, df_unempl, df_mortality)

## Create initial analysis dataset for replication 
df_09_17 <- df %>%
  mutate(aca_expanders = ifelse(state_fips %in% aca, 1, 0)) %>% ## Assign treatment
  mutate(treated = ifelse(aca_expanders == 1 & year >= 2014, 1, 0)) %>% ## Assign treatment
  filter(!(state_fips %in% drop)) %>% ## Drop states with later expansions 
  filter(year >= 2009) %>% filter(year <= 2017)

## Set of missing-data-related sample restrictions 
df_09_17 <- df_09_17 %>%
  filter(!is.na(deaths)) %>% ## Missing mortality data (9, one VA county 51917)
  filter(!is.na(unemp))  ## Missing Unemployment (227)

## Grab sample 
sample_09_17 <- df_09_17 %>% 
  group_by(state_abb, county_fips, county ) %>% 
  summarise(n =n()) %>% 
  group_by()

## Merge back on to drop unbalanced (None dropped)
df_09_17 <- df_09_17 %>% 
  left_join(sample_09_17, join_by(state_abb, county_fips, county )) %>%
  filter(n == 9) %>%
  select(-n)

#### FUNCTIONS 



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
      ungroup() %>%
      as.numeric()
    
    return(ate)
    
  }
  
}


## DEFINE SDID FUNCTION 
f_sdid <- function(data, 
                   var_state, 
                   var_period, 
                   var_y, 
                   var_did,
                   weighted = 0, 
                   weight_var = FALSE, 
                   covar_list = FALSE, # MAKE EMPTY LIST
                   calc_se = 0, # 0 = no SEs, 1 = calculate SEs
                   n_boot = 100
) {
  
  
  ## COVARIATES 
  if (covar_list != FALSE) {
    
  }
  else {
    
    # print("No covariates")
    
  }
  
  ## OPTION 1: No weights   
  if (weighted == 0) { 
    
    ## Display weights 
    print("Run without weights")
    
    ## Basic set-up 
    setup <- data %>%
      mutate(state = !! sym(var_state), 
             period = !! sym(var_period),
             y = !! sym(var_y),
             did = !! sym(var_did)) %>%
      dplyr::select(state, period, y, did) %>%
      as.data.frame() %>%
      synthdid::panel.matrices(
        unit = 1,
        time = 2,
        outcome = 3,
        treatment = 4
      )
    
    ## Run SDID 
    est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
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
      mutate(weight = !! sym(weight_var), 
             state = !! sym(var_state), 
             period = !! sym(var_period),
             y = !! sym(var_y),
             did = !! sym(var_did)) %>%
      select(collapse_var, period, weight, state, period, did, treated) %>%
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
      mutate(weight = !! sym(weight_var), 
             state = !! sym(var_state), 
             period = !! sym(var_period),
             y = !! sym(var_y),
             did = !! sym(var_did)) %>%
      select(collapse_var, period, weight, state, period, did, treated) %>%
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


## DEFINE DID FUNCTION 
run_did_function <- function(data, 
                             weights_option = FALSE,
                             weights_var = NA, 
                             controls_option = TRUE,
                             n = 1,
                             model_list_name = models){
  
  ## Set up temporary data file 
  tmp <- data 
  
  ## SET UP WEIGHTS 
  
  ## IF true, let weights = weights_var
  if(weights_option == TRUE) {
    tmp <- tmp %>% mutate(weight = pop_20_64)
  }
  
  ### ELSE false, let weights = 1 
  if(weights_option == FALSE) {
    tmp <- tmp %>% mutate(weight = 1)
  }
  
  ## SET UP CONTROLS 
  if(controls_option == TRUE) {
    f <- paste0("crude_rate ~ treated + pct_white + pct_55_64 + log_20_64 + log_35_44 + log_f_20_64 + unemp | county_fips + year")
  }
  
  if(controls_option == FALSE) {
    f <- paste0("crude_rate ~ treated  | county_fips + year")
  }
  
  ## SET UP FORMULA 
  f <- as.formula(f)
  
  ## RUN MODEL
  m <- feols(f,
             se = "cluster", 
             cluster = "state_fips", 
             weights = ~ weight,
             data = tmp)
  
  return(m)
  
}

### PRODUCE ANALYSIS DATASETS 
df_sdid <- df_09_17 %>%
  select(county_fips, year, crude_rate, treated) 


### RUN MODELS 

tmp <- f_sdid(df_sdid, "county_fips", "year", "crude_rate", "treated", weighted = 0)

## Create dataset with required variables 
df_sdid_data <- df_09_17 %>%
  select(county_fips, year, crude_rate, treated) 

## Convert to data frame 
df_sdid_data <- as.data.frame(df_sdid_data)

## Run panel matrix setup command 
setup <- panel.matrices(df_sdid_data)

## Get matrix of covariates 
tmp <- df_09_17 %>%
  select(county_fips, year, pct_white, pct_55_64, log_20_64, log_35_44, log_f_20_64, unemp) %>%
  array(dim = dim = c(dim(setup$Y))) 

X <- array(data = tmp, dim = c(dim(setup$Y)))

est <- synthdid_estimate(setup$Y, setup$N0, setup$T0, X)
plot(est)


## CREAT MODEL LIST 
models <- list()

models[[paste("nw", "nc", sep = ".")]] <- 
  run_did_function(df_09_17, 
                 weights_option = FALSE, 
                 controls_option = FALSE, 
                 model_list_name = models)

models[[paste("w", "nc", sep = ".")]] <- 
  run_did_function(df_09_17, 
                 weights_option = TRUE, 
                 weights_var = pop_20_64, 
                 controls_option = FALSE, 
                 model_list_name = models)

models[[paste("nw", "c", sep = ".")]] <- 
  run_did_function(df_09_17, 
                   var_state = county_fips, 
                   var_period = year, 
                   var_y = crude_rate, 
                   var_did = treated,
                   weights_option = FALSE, 
                   controls_option = TRUE, 
                   model_list_name = models)

models[[paste("w", "c", sep = ".")]] <- 
  run_did_function(df_09_17, 
                   weights_option = TRUE, 
                   weights_var = pop_20_64, 
                   controls_option = TRUE, 
                   model_list_name = models)


## TEST NAMING 
v = c("ATE")
names(v) = c("treatedTRUE")

modelsummary(models, output = "gt", 
             stars = c('*' = .1, '**' = .05, '***' = .01), 
             coef_map = v, 
             coef_rename = c('nw.nc' = 'Test'))  %>% 
  # column labels
  tab_spanner(label = 'Without Controls', columns = 2:3) %>%
  tab_spanner(label = 'With Controls', columns = 4:5)


msummary(models, stars = c('*' = .1, '**' = .05, '***' = .01), 
         coef_map = v, 
         shape = "rbind")
                                                                                                           model2), "Panel B" = list(model3, model4)))

names(v)





msummary(did_w_nc, stars = c('*' = .1, '**' = .05, '***' = .01))

# feols clusters by the first
# fixed effect by default, no adjustment necessary
did_nw_nc <- feols(crude_rate ~ treated | county_fips + year,
              se = "cluster", cluster = "state_fips", 
              data = df_09_17)

models[["No weights, no controls"]] <- did_nw_nc




# Interact quarter with being in the treated group using
# the fixest i() function, which also lets us specify
# a reference period (using the numeric version of Quarter)
clfe <- feols(crude_rate ~ i(year, aca_expanders, ref = 2013) | county_fips + year,
              se = "cluster", cluster = "state_fips", 
              weights = ~ pop_20_64,
              data = df_09_17)

# And use coefplot() for a graph of effects
coefplot(clfe)

clfe <- feols(crude_rate ~ i(year, aca_expanders, ref = 2013) | county_fips + year,
              se = "cluster", cluster = "state_fips", 
              data = df_09_17)

# And use coefplot() for a graph of effects
coefplot(clfe)


boot_strat_newID <- function(df_ids, # dataset with IDs and state of every station
                             df_full, # full dataset
                             strat_n, # named vector giving # of unique stations in every state (names must be state names)
                             strat_var, # name of variable you're stratifying on
                             seed = 1234){ # need seed argument so you can set inside bootstrap loop to fix draws that happen in other functions/code
  set.seed(seed)
  ids <- foreach(i = 1:length(strat_n), .combine = c) %do% {
    return(sample(df_ids[, "id"][df_ids[, strat_var] == names(strat_n)[i]], 
                  strat_n[i], replace = T))}
  df_full %<>%
    left_join(data.frame(id = ids,
                         boot_id = 1:length(ids)),
              by = "id",
              multiple = "all") %>%
    filter(!is.na(boot_id)) 
  df_full[, !ls(df_full) %in% "id"] %>% 
    rename(id = boot_id) %>%
    return
}










# helper functions for constructing factors from the sample names and classifying states based on boostrapped fits

# clean up the sample name into a nice string, and order them for plotting 
add_sample_factor <- function(x, main_spec){
  mutate(x, 
         sample_factor = ifelse(grepl("regional_", sample_name), "Regional breaks", sample_name) %>% 
           gsub("drop", "Drop ", .) %>% 
           {paste0(., ifelse(grepl("_y", .), " years", ""))} %>% 
           gsub("_y", " obs, ", .) %>% 
           gsub("^n", "", .) %>% 
           # for the main spec, label it as such 
           {ifelse(sample_name == main_spec, paste0("Main (", ., ")"), .)} %>% 
           # for a couple different cases, manually specify the clean string for the specification name
           {case_match(., 
                       "all" ~ "All obs", 
                       "discontinuous" ~ "Discontinuous", 
                       "pos_anom" ~ "Positive smoke\nPM2.5 anomalies", 
                       .default = .)} %>% 
           # now order them for plotting 
           factor(levels = c("Main (50 obs, 15 years)", 
                             "All obs", 
                             "50 obs, 5 years", 
                             "50 obs, 10 years", 
                             "100 obs, 15 years", 
                             "Drop 2020", 
                             "Drop 2021", 
                             "Drop 2022", 
                             "Regional breaks", 
                             "Discontinuous", 
                             "Positive smoke\nPM2.5 anomalies"), 
                  ordered = T))
}

# using the distribution of coefficient estimates, calculate quantiles and classify as smoke-influenced or not 
classify_boot <- function(grouped_df, full = T){
  # for all of the non-grouping columns, calculate the median, alpha/2 and 1 - alpha/2 quantiles 
  out <- summarise(grouped_df, 
                   across(everything(), 
                          list(lower =~ quantile(.x, q_alpha/2), 
                               mid =~ quantile(.x, 0.5),
                               upper =~ quantile(.x, 1 - q_alpha/2), 
                               # if both alpha/2 and 1 - alpha/2 quantiles have the same sign, its siggy 
                               sig =~ sign(quantile(.x, 1 - q_alpha/2)) == sign(quantile(.x, q_alpha/2)))),
                   n_boot = n(),
                   .groups = "drop") %>% 
    # define it as smoke influenced if the upper limit of b2nonsmoke - total < 0
    mutate(smoke_influenced = b2nonsmokeMb2total_upper < 0, 
           smoke_influenced_text = ifelse(smoke_influenced,
                                          "smoke-influenced",
                                          "no smoke influence detected"))
  if(full){ # if specified, do the full classification (for average PM2.5 and not extreme daily PM2.5)
    out %<>% mutate(early_decline = b1total_upper < 0, 
                    # classify the change in the total PM trend
                    change_group = case_when(early_decline & b2total_lower > 0 ~ "reversal",
                                             early_decline & b2totalMb1total_upper < 0 ~ "acceleration",
                                             early_decline & b2totalMb1total_lower > 0 ~ "stagnation",
                                             early_decline & !b2totalMb1total_sig ~ "non-sig change",
                                             !early_decline ~ "no sig early decline",
                                             T ~ "something wrong"), 
                    # classify the difference between the non-smoke and total PM trends in the post period 
                    smoke_group = case_when(!smoke_influenced ~ "no smoke influence detected", 
                                            smoke_influenced & !early_decline ~ "smoke-influenced,\nno early decline",
                                            # isolate smoke-casued reversals first 
                                            smoke_influenced & 
                                              b2total_lower > 0 & b2nonsmoke_upper < 0 ~ "smoke-caused reversal",
                                            # remaining reversals are just smoke-influenced
                                            smoke_influenced & 
                                              b2total_lower > 0 ~ "smoke-influenced reversal", 
                                            # smoke-influenced and stagnating 
                                            smoke_influenced & 
                                              b2totalMb1total_lower > 0 ~ "smoke-influenced stagnation", 
                                            # what's left could either not be stagnating (would have improved faster/accelerated without smoke)
                                            # or had no early decline in total PM
                                            smoke_influenced & b2totalMb1total_upper < 0 ~ "smoke-influenced acceleration",
                                            smoke_influenced & !b2totalMb1total_sig ~ "smoke-influenced,\nno sig. trend change",
                                            T ~ "something wrong")) %>% 
      # now make ordered factor for both classifications, again for conventient plotting
      mutate(smoke_group_factor = factor(smoke_group,
                                         levels = c("smoke-caused reversal",
                                                    "smoke-influenced reversal",
                                                    "smoke-influenced stagnation",
                                                    "smoke-influenced,\nno early decline",
                                                    "smoke-influenced,\nno sig. trend change",
                                                    "smoke-influenced acceleration",
                                                    "no smoke influence detected"),
                                         ordered = TRUE),
             change_group = gsub("sig", "sig.", change_group),
             change_group_factor = factor(change_group,
                                          levels = c("reversal",
                                                     "stagnation",
                                                     "acceleration",
                                                     "non-sig. change",
                                                     "no sig. early decline"),
                                          ordered = TRUE))
  }
  # also add state abbreviation to the df using a state dictionary 
  out %>% 
    rename(state_name = state) %>% 
    mutate(state_name = case_match(state_name,
                                   "US" ~ "Contiguous US", 
                                   .default = state_name)) %>% 
    left_join(state_dict) %>% 
    mutate(state_abbr = case_match(state_name, 
                                   "Contiguous US" ~ "US", 
                                   .default = state_abbr)) %>% 
    return
} 










classify_trends <- function(df){
  df %>%
    mutate(early_decline = (b1total_est < 0 & b1total_pval < 0.05), 
           smoke_influenced = (b2nonsmokeMb2total_est < 0 & b2nonsmokeMb2total_pval < 0.05),
           change_group = case_when(early_decline & b2total_est > 0 & b2total_pval < 0.05 ~ "reversal",
                                    early_decline & b2totalMb1total_est < 0 & b2totalMb1total_pval < 0.05 ~ "acceleration",
                                    early_decline & b2totalMb1total_est > 0 & b2totalMb1total_pval < 0.05 ~ "stagnation",
                                    early_decline & b2totalMb1total_pval >= 0.05 ~ "non-sig change",
                                    !early_decline ~ "no sig early decline",
                                    T ~ "something wrong"), 
           smoke_group = case_when(!smoke_influenced ~ "no smoke influence detected", 
                                   # isolate smoke-casued reversals first 
                                   early_decline & smoke_influenced & 
                                     b2total_est > 0 & b2total_pval < 0.05 & 
                                     b2nonsmoke_est < 0 & b2nonsmoke_pval < 0.05 ~ "smoke-caused reversal",
                                   # remaining reversals are just smoke-influenced
                                   early_decline & smoke_influenced & 
                                     b2total_est > 0 & b2total_pval < 0.05 ~ "smoke-influenced reversal", 
                                   # smoke-influenced and stagnating 
                                   early_decline & smoke_influenced & 
                                     b2totalMb1total_est > 0 & b2totalMb1total_pval < 0.05 ~ "smoke-influenced stagnation", 
                                   # what's left could either not be stagnating (would have improved faster/accelerated without smoke)
                                   # or had no early decline in total PM
                                   !early_decline & smoke_influenced ~ "smoke-influenced, no early decline",
                                   smoke_influenced & b2totalMb1total_pval >= 0.05 ~ "smoke-influenced, no signficant trend change",
                                   T ~ "something wrong")) %>% 
    return
}










# function for cleaning CONUS small multiples----
clean_geo_facet <- function(fig){
  # add x axes to bottom panel in each column
  fig %<>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-1-8"),
                    t = 11, l = 5, b = 11, r = 5,
                    name = "US-axis-b") %>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-1-8"),
                    t = 32, l = 5, b = 32, r = 5,
                    name = "CA-axis-b") %>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-2-8"),
                    t = 37, l = 9, b = 37, r = 9,
                    name = "AZ-axis-b") %>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-3-8"),
                    t = 37, l = 13, b = 37, r = 13,
                    name = "NM-axis-b") %>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-5-8"),
                    t = 42, l = 21, b = 42, r = 21,
                    name = "LA-axis-b") %>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-6-8"),
                    t = 42, l = 25, b = 42, r = 25,
                    name = "MS-axis-b") %>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-7-8"),
                    t = 42, l = 29, b = 42, r = 29,
                    name = "AL-axis-b") %>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-9-8"),
                    t = 32, l = 37, b = 32, r = 37,
                    name = "MD-axis-b") %>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-10-8"),
                    t = 32, l = 41, b = 32, r = 41,
                    name = "DE-axis-b") %>%
    gtable_add_grob(gtable_filter(fig, pattern = "axis-b-9-8"),
                    t = 27, l = 45, b = 27, r = 45,
                    name = "RI-axis-b")
  
  # remove extra axes
  fig %<>% gtable_remove_grobs(c(paste0("axis-l-1-", 2:10),
                                 paste0("axis-l-2-", 1:5),
                                 # paste0("strip-t-2-", 1:5),
                                 paste0("axis-l-2-", 7:9),
                                 paste0("axis-l-3-", c(8, 11)),
                                 paste0("axis-l-6-", c(1, 9:11)),
                                 paste0("axis-l-7-", c(1:3, 9:11)),
                                 paste0("axis-l-8-", c(1:3,5:7,9:11)),
                                 paste0("axis-b-", c(1:3, 5:7, 9:11), "-8"),
                                 # wildly nonsensical panel and strip naming scheme
                                 "strip-t-1-2", "strip-t-2-2", "strip-t-2-1",
                                 "panel-2-1", "panel-9-1", "panel-10-1"))
  # make the US panel bigger
  fig$layout %<>% mutate(b = ifelse(name %in% c("panel-1-1", "axis-l-1-1"),
                                    12.9, b),
                         b = ifelse(name %in% c("US-axis-b"),
                                    13.9, b),
                         t = ifelse(name %in% c("US-axis-b"),
                                    12.9, t),
                         r = ifelse(name %in% c("panel-1-1", "strip-t-1-1", "US-axis-b"), # "axis-b-1-1",
                                    8.9, r))
  
  # add y axis label to US panel as well
  fig %<>% gtable_add_grob(gtable_filter(fig, pattern = "ylab-l"),
                           t = 8, l = 3, b = 12.9, r = 3,
                           name = "ylab-l-US")
  return(fig)
}

# set up state grid
my_grid <- geofacet::us_state_grid3
# move ME up, VT and NH one right, FL one left
my_grid %<>% 
  mutate(row = ifelse(code == "ME", row - 1, row), 
         col = ifelse(code %in% c("VT", "NH"), col + 1, col), 
         col = ifelse(code == "FL", col - 1, col)) %>% 
  # remove AK and HI
  filter(code %in% c("AK", "HI") == FALSE) %>% 
  # add US to top left
  rbind(data.frame(row = 1, col = 1, code = "US", name = "Contiguous US"))










fit_break_trend <- function(station_year_df, 
                            fit_family_fn = feglm,
                            PM_measure = c("mean", "extremeN", "extremePct"), #
                            break_alpha = 0.1, 
                            break_var = "year", 
                            tol = 1e-7, # convergence tolerance for breakpoint algorithm
                            max_it = 10, 
                            extra_output = F,# keep this off by default, will break estimation code
                            corner_push = 0.95, # scalar that governs how far to push algorithm off of corner solution, lower = bigger adjustment
                            max_restart = 3, # maximum number of times to restart algorithm
                            fe_fml = NULL, # character, what FE you want to use, should have any operators like + or ^ built in
                            break_init = NULL, # should be a vector of length > 1 to build a uniform distribution with. Recommended that this comes from the output of the smart_initialize function
                            offset = NULL, 
                            glm_it = 15, 
                            ...){ 
  # condition checks
  if(is.null(fe_fml)){fit_fml <- as.formula(paste0("value ~ years_pre + years_post"))}
  if(!is.null(fe_fml)){fit_fml <- as.formula(paste0("value ~ years_pre + years_post | ", fe_fml))}
  if(PM_measure == "extremeN"){family <- "poisson"}
  if(PM_measure == "extremePct"){family <- "binomial"}
  if(PM_measure == "mean"){family <- "gaussian"}
  
  ## reshape data
  station_year_PM <- station_year_df %>% # pivot the station-year df long, so its station-year-PM type (i.e. non-smoke vs total)
    pivot_longer(contains("PM"),
                 names_sep = "_",
                 names_to = c("pm_type", ".value")) %>% 
    dplyr::rename(value = all_of(PM_measure)) %>% 
    filter(!is.na(value)) 
  
  
  ## estimation
  break_mod <- fixest_segmented(station_year_PM, tol = tol, 
                                fit_family_fn = fit_family_fn, 
                                Z_var = break_var, alpha = break_alpha, 
                                max_it = max_it, fe_fml = fe_fml, 
                                fc = corner_push, max_restart = max_restart, 
                                psi0 = break_init, offset = offset, 
                                family = family)
  break_year <- break_mod$psi
  
  # using estimated break, fit state specific trends 
  est_df <- station_year_PM %>% 
    mutate(period = (year > break_year),
           years_post = (year > break_year)*(year - break_year), 
           years_pre = (year <= break_year)*(year - break_year)) 
  
  df_out <- est_df %>% 
    {rbind(.,
           mutate(., state = "US"))} %>%
    nest_by(state) %>%
    mutate(sd = data %>% 
             mutate(period = ifelse(period, "post", "pre")) %>% 
             summarise(sd = sd(value, na.rm = T), 
                       .by = c(period, pm_type)) %>% 
             pivot_wider(values_from = sd, names_from = c(period, pm_type), 
                         names_prefix = "sd_") %>% 
             unlist %>% list,
           mean = data %>% 
             mutate(period = ifelse(period, "post", "pre")) %>% 
             summarise(mean = mean(value, na.rm = T), 
                       .by = c(period, pm_type)) %>% 
             pivot_wider(values_from = mean, names_from = c(period, pm_type), 
                         names_prefix = "mean_") %>% 
             unlist %>% list,
           mod_total = list(try(fit_family_fn(fit_fml,
                                              data = data %>% filter(pm_type == "totalPM"),
                                              glm.iter = glm_it, 
                                              fixef.rm = fixef_rm, 
                                              family = family,
                                              offset = offset, ...),
                                silent = T)), # always set threads to 1 here since you're going to do be doing this in parallel
           mod_nonsmoke = list(try(fit_family_fn(fit_fml,
                                                 data = data %>% filter(pm_type == "nonsmokePM"),
                                                 glm.iter = glm_it, 
                                                 fixef.rm = fixef_rm, 
                                                 family = family, 
                                                 offset = offset, ...),
                                   silent = T)),
           nobs_total = if(class(mod_total) == "fixest"){nobs(mod_total)} else{NA},
           nobs_nonsmoke = if(class(mod_nonsmoke) == "fixest"){nobs(mod_nonsmoke)} else{NA},
           it_total = if(class(mod_total) == "fixest"){mod_total$iterations} else{NA},
           conv_total = if(class(mod_total) == "fixest"){mod_total$convStatus} else{NA},
           it_nonsmoke = if(class(mod_nonsmoke) == "fixest"){mod_nonsmoke$iterations} else{NA}, 
           conv_nonsmoke = if(class(mod_nonsmoke) == "fixest"){mod_nonsmoke$convStatus} else{NA},
           b1total = if(class(mod_total) == "fixest"){coef(mod_total) %>% magrittr::extract("years_pre")} else{NA}, 
           b2total = if(class(mod_total) == "fixest"){coef(mod_total) %>% magrittr::extract("years_post")} else{NA}, 
           b1nonsmoke = if(class(mod_nonsmoke) == "fixest"){coef(mod_nonsmoke) %>% magrittr::extract("years_pre")} else{NA}, 
           b2nonsmoke = if(class(mod_nonsmoke) == "fixest"){coef(mod_nonsmoke) %>% magrittr::extract("years_post")} else{NA},
           break_year = break_year, 
           it_algo = break_mod$iterations) %>% 
    dplyr::select(-data, -starts_with("mod")) %>% 
    unnest_wider(c(sd, mean))
  
  if(extra_output){ # set this to T only for debugging/checking results, this will break rest of the pipeline in normal code
    return(list(df_out = df_out, 
                L_list = break_mod$L_list, 
                psi_list = break_mod$psi_list, 
                eps_list = break_mod$eps_list, 
                mod_list = break_mod$mod_list,
                iterations = break_mod$iterations))}
  if(!extra_output){
    return(df_out)
  }
}










fit_break_discontinuous <- function(station_year_df, 
                                    fit_family_fn = femlm, # 
                                    break_fitstat = "wll", # r2, wr2 (within r2), cor2, ll, wll (within ll), delta_b
                                    PM_measure = NULL, # should be "mean", "extremePct" or "extremeN"
                                    break_alpha = 0.2, 
                                    breaks_by = 0.5, # default to half-year breaks 
                                    offset = NULL, 
                                    fixef_rm = "perfect", 
                                    ...){ 
  
  fitstat_register("wll", function(x) fitstat(x, "ll", simplify = T) - x$ll_fe_only , "within LL")
  
  # condition checks
  if(!PM_measure %in% c("mean", "extremeN")){stop("PM_measure can only be 'mean' or 'extremeN'")}
  if(PM_measure == "mean"){family <- "gaussian"}else{family <- "poisson"}
  if(break_fitstat %in% c("rmse", "dev")){fitstat_dir = 1}
  if(break_fitstat %in% c("r2", "wll", "ll")){fitstat_dir = -1}
  
  # define the FE based on whether its piecewise continuous or discontinuous
  fe_fml <- "id^period" 
  fit_fml <- as.formula(paste0("value ~ years_pre + years_post | ", fe_fml))
  
  # pivot the station-year df long, so its station-year-PM type (i.e. non-smoke vs total)
  
  station_year_PM <- station_year_df %>% # pivot the station-year df long, so its station-year-PM type (i.e. non-smoke vs total)
    pivot_longer(contains("PM"),
                 names_sep = "_",
                 names_to = c("pm_type", ".value")) %>% 
    dplyr::rename(value = all_of(PM_measure)) %>% 
    filter(!is.na(value)) 
  
  # possible set of years for the break to happen at 
  opt_breaks <- station_year_PM$year %>% 
    unique %>% 
    sort %>% 
    quantile(c(break_alpha, 1 - break_alpha)) %>% 
    # round lower bound up, upper bound down
    {seq(from = ceiling(.[1]), 
         to = floor(.[2]), 
         by = breaks_by)} 
  # for each possible break year, fit a US wide model
  break_year <- purrr::map_dbl(opt_breaks, function(b_year){
    fit_family_fn(fit_fml, 
                  offset = offset, 
                  family = family,
                  data = station_year_PM %>% 
                    filter(pm_type == "totalPM") %>% 
                    mutate(period = (year > b_year),
                           years_post = (year > b_year)*(year - b_year), 
                           years_pre = (year <= b_year)*(year - b_year)), ...) %>% 
      fitstat(break_fitstat, simplify = T)}) %>% 
    multiply_by(fitstat_dir) %>%
    which.min() %>% 
    magrittr::extract(opt_breaks, .) 
  
  # using estimated break year, fit models
  est_df <- station_year_PM %>% 
    mutate(period = (year > break_year),
           years_post = (year > break_year)*(year - break_year), 
           years_pre = (year <= break_year)*(year - break_year)) 
  
  df_out <- est_df %>% 
    {rbind(.,
           mutate(., state = "US"))} %>%
    nest_by(state) %>%
    mutate(sd = data %>% 
             mutate(period = ifelse(period, "post", "pre")) %>% 
             summarise(sd = sd(value, na.rm = T), 
                       .by = c(period, pm_type)) %>% 
             pivot_wider(values_from = sd, names_from = c(period, pm_type), 
                         names_prefix = "sd_") %>% 
             unlist %>% list,
           mean = data %>% 
             mutate(period = ifelse(period, "post", "pre")) %>% 
             summarise(mean = mean(value, na.rm = T), 
                       .by = c(period, pm_type)) %>% 
             pivot_wider(values_from = mean, names_from = c(period, pm_type), 
                         names_prefix = "mean_") %>% 
             unlist %>% list,
           mod_total = list(try(fit_family_fn(fit_fml,
                                              data = data %>% filter(pm_type == "totalPM"),
                                              fixef.rm = fixef_rm, 
                                              family = family,
                                              offset = offset, ...),
                                silent = T)), # always set threads to 1 here since you're going to do be doing this in parallel
           mod_nonsmoke = list(try(fit_family_fn(fit_fml,
                                                 data = data %>% filter(pm_type == "nonsmokePM"),
                                                 fixef.rm = fixef_rm, 
                                                 family = family, 
                                                 offset = offset, ...),
                                   silent = T)),
           nobs_total = if(class(mod_total) == "fixest"){nobs(mod_total)} else{NA},
           nobs_nonsmoke = if(class(mod_nonsmoke) == "fixest"){nobs(mod_nonsmoke)} else{NA},
           it_total = if(class(mod_total) == "fixest"){mod_total$iterations} else{NA},
           conv_total = if(class(mod_total) == "fixest"){mod_total$convStatus} else{NA},
           it_nonsmoke = if(class(mod_nonsmoke) == "fixest"){mod_nonsmoke$iterations} else{NA}, 
           conv_nonsmoke = if(class(mod_nonsmoke) == "fixest"){mod_nonsmoke$convStatus} else{NA},
           b1total = if(class(mod_total) == "fixest"){coef(mod_total) %>% magrittr::extract("years_pre")} else{NA}, 
           b2total = if(class(mod_total) == "fixest"){coef(mod_total) %>% magrittr::extract("years_post")} else{NA}, 
           b1nonsmoke = if(class(mod_nonsmoke) == "fixest"){coef(mod_nonsmoke) %>% magrittr::extract("years_pre")} else{NA}, 
           b2nonsmoke = if(class(mod_nonsmoke) == "fixest"){coef(mod_nonsmoke) %>% magrittr::extract("years_post")} else{NA},
           break_year = break_year) %>% 
    dplyr::select(-data, -starts_with("mod")) %>% 
    unnest_wider(c(sd, mean))
  
  return(df_out)
  
  
}










fit_break_trend <- function(station_year_df, 
                            fit_family_fn = feglm,
                            PM_measure = c("mean", "extremeN"), #
                            break_alpha = 0.2, 
                            break_var = "year", 
                            tol = 1e-8, # convergence tolerance for breakpoint algorithm
                            max_it = 5, 
                            extra_output = F,# keep this off by default, will break estimation code
                            corner_push = 0.95, # scalar that governs how far to push algorithm off of corner solution, lower = bigger adjustment
                            max_restart = 2, # maximum number of times to restart algorithm
                            fe_fml = NULL, # character, what FE you want to use, should have any operators like + or ^ built in
                            break_init = NULL, # should be a vector of length > 1 to build a uniform distribution with. Recommended that this comes from the output of the smart_initialize function
                            offset = NULL, 
                            glm_it = 15, 
                            ...){ 
  # condition checks
  if(is.null(fe_fml)){fit_fml <- as.formula(paste0("value ~ years_pre + years_post"))}
  if(!is.null(fe_fml)){fit_fml <- as.formula(paste0("value ~ years_pre + years_post | ", fe_fml))}
  if(PM_measure == "extremeN"){family <- "poisson"}
  if(PM_measure == "mean"){family <- "gaussian"}
  
  ## reshape data
  station_year_PM <- station_year_df %>% # pivot the station-year df long, so its station-year-PM type (i.e. non-smoke vs total)
    pivot_longer(contains("PM"),
                 names_sep = "_",
                 names_to = c("pm_type", ".value")) %>% 
    dplyr::rename(value = all_of(PM_measure)) %>% 
    filter(!is.na(value)) 
  
  
  ## estimation
  break_mod <- fixest_segmented(station_year_PM, tol = tol, 
                                fit_family_fn = fit_family_fn, 
                                Z_var = break_var, alpha = break_alpha, 
                                max_it = max_it, fe_fml = fe_fml, 
                                fc = corner_push, max_restart = max_restart, 
                                psi0 = break_init, offset = offset, 
                                family = family)
  break_year <- break_mod$psi
  
  # using estimated break, fit state specific trends 
  est_df <- station_year_PM %>% 
    mutate(period = (year > break_year),
           years_post = (year > break_year)*(year - break_year), 
           years_pre = (year <= break_year)*(year - break_year)) 
  
  df_out <- est_df %>% 
    {rbind(.,
           mutate(., state = "US"))} %>%
    nest_by(state) %>%
    mutate(sd = data %>% 
             mutate(period = ifelse(period, "post", "pre")) %>% 
             summarise(sd = sd(value, na.rm = T), 
                       .by = c(period, pm_type)) %>% 
             pivot_wider(values_from = sd, names_from = c(period, pm_type), 
                         names_prefix = "sd_") %>% 
             unlist %>% list,
           mean = data %>% 
             mutate(period = ifelse(period, "post", "pre")) %>% 
             summarise(mean = mean(value, na.rm = T), 
                       .by = c(period, pm_type)) %>% 
             pivot_wider(values_from = mean, names_from = c(period, pm_type), 
                         names_prefix = "mean_") %>% 
             unlist %>% list,
           mod_total = list(try(fit_family_fn(fit_fml,
                                              data = data %>% filter(pm_type == "totalPM"),
                                              glm.iter = glm_it, 
                                              fixef.rm = fixef_rm, 
                                              family = family,
                                              offset = offset, ...),
                                silent = T)), # always set threads to 1 here since you're going to do be doing this in parallel
           mod_nonsmoke = list(try(fit_family_fn(fit_fml,
                                                 data = data %>% filter(pm_type == "nonsmokePM"),
                                                 glm.iter = glm_it, 
                                                 fixef.rm = fixef_rm, 
                                                 family = family, 
                                                 offset = offset, ...),
                                   silent = T)),
           nobs_total = if(class(mod_total) == "fixest"){nobs(mod_total)} else{NA},
           nobs_nonsmoke = if(class(mod_nonsmoke) == "fixest"){nobs(mod_nonsmoke)} else{NA},
           it_total = if(class(mod_total) == "fixest"){mod_total$iterations} else{NA},
           conv_total = if(class(mod_total) == "fixest"){mod_total$convStatus} else{NA},
           it_nonsmoke = if(class(mod_nonsmoke) == "fixest"){mod_nonsmoke$iterations} else{NA}, 
           conv_nonsmoke = if(class(mod_nonsmoke) == "fixest"){mod_nonsmoke$convStatus} else{NA},
           b1total = if(class(mod_total) == "fixest"){coef(mod_total) %>% magrittr::extract("years_pre")} else{NA}, 
           b2total = if(class(mod_total) == "fixest"){coef(mod_total) %>% magrittr::extract("years_post")} else{NA}, 
           b1nonsmoke = if(class(mod_nonsmoke) == "fixest"){coef(mod_nonsmoke) %>% magrittr::extract("years_pre")} else{NA}, 
           b2nonsmoke = if(class(mod_nonsmoke) == "fixest"){coef(mod_nonsmoke) %>% magrittr::extract("years_post")} else{NA},
           break_year = break_year, 
           it_algo = break_mod$iterations) %>% 
    dplyr::select(-data, -starts_with("mod")) %>% 
    unnest_wider(c(sd, mean))
  
  if(extra_output){ # set this to T only for debugging/checking results, this will break rest of the pipeline in normal code
    return(list(df_out = df_out, 
                L_list = break_mod$L_list, 
                psi_list = break_mod$psi_list, 
                eps_list = break_mod$eps_list, 
                mod_list = break_mod$mod_list,
                iterations = break_mod$iterations))}
  if(!extra_output){
    return(df_out)
  }
}










# search function to be optimized
search_min <- function(h, psi, psi_old, dt, 
                       fit_family_fn, 
                       family = NULL, 
                       offset = NULL, 
                       ...) {
  
  
  psi_ok <- psi*h + psi_old*(1 - h)
  obj1 <- fit_family_fn(value ~ Z + u,
                        data = dt %>%
                          mutate(u = (Z > psi_ok)*(Z - psi_ok)), 
                        family = family, 
                        offset = offset, 
                        ...)
  if("deviance" %in% names(obj1)){
    L1 <- obj1$deviance
  } else {L1 <- sum(obj1$residuals^2)}
  return(L1)
}

adj_psi <- function(psii, lim) {pmin(pmax(lim[1],psii),lim[2])} 

# function to push algorithm away from corner solutions
psi_corner <- function(psii, lim, fc){
  if(psii <= lim[1]){
    return(lim[1] + diff(lim)*(1-fc))}
  if(psii >= lim[2]){
    return(lim[2] - diff(lim)*(1-fc))}
  else{return(psii)}    # if psi <= lower lim, multiply the range by 1/fc and add. if psi >= upper lim, multiply the range by fc and subtract, 
  # if lower lim < psi < upper lim, do nothing
}

# fitting function, with some updates 
fixest_segmented <- function(dt = station_year_PM, # data frame 
                             tol = 1e-8, # convergence tolerance
                             fit_family_fn = feglm, 
                             Z_var = "year", # breakpoint variable, character
                             alpha = 0.2, # restrict breakpoint search to (alpha, 1-alpha) quantiles
                             psi0 = NULL, # initialization value, should be a sequence not a single value
                             max_it = 5, # max iterations
                             max_restart = 2, # how many restarts
                             fc = 0.95, # inflation/deflation factor for when algorithm approaches edge of search space
                             fe_fml = NULL, # character name of formula, should include operators like + or ^ since it's just a straight paste, no collapse
                             dep_var = "value", # name of variable
                             return_all = F, # set this to T for more detailed output, but will increase computation time and may break piped code
                             offset = NULL, # used for poisson regressions
                             fixef_rm = "perfect", # how to deal with observations where the unit (station) has all zeroes/constant outcomes. "perfect" is fixest default
                             glm_iter = 15, 
                             family = "gaussian",
                             ...){
  # check conditions
  if(!is.null(fe_fml)){est_eq <- as.formula(paste(dep_var, "~ Z + u + v |", fe_fml))}
  if(is.null(fe_fml)){est_eq <- as.formula(paste(dep_var, "~ Z + u + v"))} # estimation equation
  if(length(psi0) == 1){stop("psi0 must be a sequence of length 2 or greater.")} # if you feed it only one value, throw an error
  
  # declare globals 
  dt %<>% rename(Z = all_of(Z_var))
  Z_range <- range(dt$Z) %>% unname
  limZ <- quantile(Z_range, probs = c(alpha, 1 - alpha)) %>% unname
  if(!is.null(psi0)){psi_new = runif(1, min(psi0), max(psi0))}else
  {psi_new = runif(1, limZ[1], limZ[2])}
  psi_list <- psi_new
  L_list <- c(NA)
  eps_list <- c(NA)
  corner_list <- c(NA)
  mod_list <- list()
  coef_list <- list()
  epsilon <- Inf
  it <- 0
  total_it <- 0
  restart <- 0
  
  # initialize
  while(abs(epsilon) > tol){
    psi <- psi_new  
    mod_fit <- fit_family_fn(est_eq,
                             data = dt %>%
                               mutate(u = (Z > psi)*(Z - psi),
                                      v = -(Z > psi)), 
                             offset = offset,
                             fixef.rm = fixef_rm, 
                             glm.iter = glm_it, 
                             family = family, ...)
    mod_coef <- mod_fit$coefficients  
    psi_new <- psi + mod_coef[["v"]]/mod_coef[["u"]]
    psi_new <- adj_psi(psi_new, Z_range) 
    a <- optimize(search_min, c(0, 1), psi = psi_new, psi_old = psi, 
                  dt = dt, fit_family_fn = fit_family_fn, 
                  family = family, offset = offset)
    use_k <- a$minimum
    L1 <- a$objective
    psi_new <- psi_new*use_k + psi*(1 - use_k)
    psi_new <- adj_psi(psi_new, limZ) 
    epsilon <- abs(psi_new - psi)/psi
    psi_in_corner <- any(psi_new == limZ)
    # if we're getting a corner solution, push in by a small amount 
    if(psi_in_corner){psi_new <- psi_corner(psi_new, limZ, fc)}
    
    coef_list <- rbind(coef_list, mod_coef)
    mod_list[[length(mod_list) + 1]] <- mod_fit
    psi_list <- c(psi_list, psi_new)
    L_list <- c(L_list, L1)
    eps_list <- c(eps_list, epsilon)
    corner_list <- c(corner_list, psi_in_corner)
    
    it <- it + 1
    total_it <- total_it + 1
    mod_final <- mod_list[length(mod_list)][[1]] # store current model in case it converges without restarts
    if (it >= max_it){
      # if you haven't met max restarts, restart everything
      if(restart < max_restart){
        restart <- restart + 1
        if(is.null(psi0)){psi_new <- runif(1, limZ[1], limZ[2])}
        if(!is.null(psi0)){psi_new <- runif(1, min(psi0), max(psi0))}
        epsilon <- Inf
        it <- 0
        print(paste0("restart ", restart, " after reaching max iterations, now starting with ", psi_new))
      } else{
        psi_new <- psi_list[which.min(L_list)]
        mod_final <- mod_list[which.min(L_list)][[1]]
        warning(paste0("Reached max iterations and max restarts without reaching tolerance.", 
                       "Selecting final psi of ", psi_new, " based on deviance value."))
        break
      }
    }
  }
  
  if(restart > 0){psi_new <- psi_list[which.min(L_list)]} # if non-zero restarts, select breakpoint based on deviance
  
  corner_solution <- ifelse(any(psi_new == limZ), T, F)
  
  if(return_all){
    return(list(psi = psi_new, 
                psi_list = psi_list, 
                L_list = L_list, 
                epsilon_list = eps_list, 
                iterations = total_it, 
                restart = restart,
                coef_list = coef_list,
                corner_list = corner_list, 
                mod_final = mod_final, 
                corner_solution = corner_solution, 
                lower_lim = limZ[1], # spit out lower and upper lims because you want to be able to figure out how close estimated psi is to corner, very close indicates that model probably thinks there is no break (at least not within the specified range)
                upper_lim = limZ[2]))}
  if(!return_all){ # return only the objects you need to speed up parallelization
    return(list(psi = psi_new, 
                iterations = total_it, 
                corner_solution = corner_solution, 
                lower_lim = limZ[1], 
                upper_lim = limZ[2]))
  }
}










smart_initialize <- function(station_year_df, 
                             fit_family_fn = feglm, 
                             break_fitstat = "r2", # r2, wr2 (within r2), cor2, ll, wll (within ll), delta_b
                             PM_measure = NULL, # should be "mean", "extremePct" or "extremeN"
                             continuous_at_break,
                             fitstat_dir = -1, # if you want to maximize, just set dir = -1
                             break_alpha = 0.2, 
                             breaks_by = 0.5, # default to half-year breaks 
                             interval_width = 1, # in index units -- how far away do you want to move from the optimal break candidate when building the interval
                             return_bounds = T, # default to T, is set to F so that function can be used to estimate discontinuous results
                             fixef_rm = "perfect", 
                             family = NULL, 
                             glm_it = 15, 
                             offset = NULL, 
                             ...){ # other args to fit_family_fn
  
  fitstat_register("dev", function(x) x$deviance)
  
  # define the FE based on whether its piecewise continuous or discontinuous
  fe_fml <- if(continuous_at_break){"id"} else{"id^period"} 
  
  # crude check to see if fitstat and its direction make sense 
  if(break_fitstat %in% c("rmse", "dev") & fitstat_dir == -1){
    warning(paste0("Do you really want to maximize the ", break_fitstat, "?"))
  } else if(break_fitstat %in% c("r2", "wr2", "cor2", "ll", "wll", "delta_b") & fitstat_dir == 1){
    warning(paste0("Do you really want to minimize the ", break_fitstat, "?"))
  }
  # pivot the station-year df long, so its station-year-PM type (i.e. non-smoke vs total)
  station_year_PM <- station_year_df %>% 
    pivot_longer(contains("PM"),
                 names_sep = "_",
                 names_to = c("pm_type", ".value")) %>% 
    rename(value = all_of(PM_measure)) %>% 
    filter(!is.na(value))
  
  # possible set of years for the break to happen at 
  opt_breaks <- station_year_PM$year %>% 
    unique %>% 
    sort %>% 
    quantile(c(break_alpha, 1 - break_alpha)) %>% 
    # round lower bound up, upper bound down
    {seq(from = ceiling(.[1]), 
         to = floor(.[2]), 
         by = breaks_by)} 
  # for each possible break year, fit a US wide model
  break_cand <- purrr::map_dbl(opt_breaks, function(b_year){
    fit_family_fn(as.formula(
      paste0("value ~ years_pre + years_post | ", fe_fml)),
      fixef.rm = fixef_rm, 
      glm.iter = glm_it, 
      family = family,
      offset = offset,
      ...,
      # paste0("value ~ years_pre + years_post | id^period")), #fe_fml)),
      data = station_year_PM %>% 
        # fit the break year using total PM only
        filter(pm_type == "totalPM") %>% 
        # define years pre/post based on proposed break year
        mutate(period = (year >= b_year),
               years_post = (year >= b_year)*(year - b_year), 
               years_pre = (year < b_year)*(year - b_year))) %>% 
      fitstat(break_fitstat, simplify = T)}) %>% 
    multiply_by(fitstat_dir)
  
  break_int <- c(which.min(break_cand) - interval_width, 
                 which.min(break_cand) + interval_width)
  break_int <- pmax(pmin(break_int, length(opt_breaks)), 1)  
  
  if(return_bounds){return(opt_breaks[break_int])}
  if(!return_bound){return(break_cand)}
  
  
}

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

fit_pooled <- function(station_year_df_wide, 
                       PM_measure, 
                       break_year, 
                       min_year_obs = 1, 
                       min_station_years = 1, 
                       y_start = -Inf, 
                       y_end = Inf,
                       y_drop = NULL, 
                       fe_regex, 
                       duplicated_break = TRUE){
  station_year_df_wide %>% 
    filter(year >= y_start & year <= y_end) %>% 
    filter(year %in% y_drop == FALSE) %>%
    mutate(period = (year > break_year),
           years_post = (year > break_year)*(year - break_year), 
           years_pre = (year <= break_year)*(year - break_year)) %>%
    filter(totalPM_n > min_year_obs) %>% 
    group_by(id) %>% 
    mutate(n_year = n()) %>% 
    ungroup %>% 
    filter(n_year >= min_station_years) %>% 
    mutate(duplicated = FALSE) %>% 
    {rbind(.,
           filter(., year == break_year) %>%
             mutate(period = TRUE,
                    duplicated = TRUE))} %>%
    filter(duplicated %in% c(FALSE, duplicated_break)) %>%
    {rbind(.,
           mutate(., state_abbr = "US"))} %>%
    pivot_longer(contains("PM"),
                 names_sep = "_",
                 names_to = c("pm_type", ".value")) %>%
    unite("fe", matches(fe_regex), remove = FALSE) %>%
    rename(value = PM_measure) %>%
    filter(!is.na(value)) %>%
    nest_by(state_abbr) %>%
    mutate(mod = list(feols(value ~ pm_type:years_pre + pm_type:years_post | fe,
                            data = data)),
           nobs = nobs(mod),
           vcv = list(vcov(mod, attr = T)),
           t_df = attr(vcv, "df.t"),
           coef = list(coef(mod)),
           b1total_inds = list((names(coef(mod)) %in% c("pm_typetotalPM:years_pre"))*1),
           b2total_inds = list((names(coef(mod)) %in% c("pm_typetotalPM:years_post"))*1),
           b1nonsmoke_inds = list((names(coef(mod)) %in% c("pm_typenonsmokePM:years_pre"))*1),
           b2nonsmoke_inds = list((names(coef(mod)) %in% c("pm_typenonsmokePM:years_post"))*1),
           b2totalMb1total_inds = list(b2total_inds - b1total_inds),
           b2nonsmokeMb2total_inds = list(b2nonsmoke_inds - b2total_inds)) %>%
    mutate(across(ends_with("inds"),
                  list(est =~ t(.x)%*%coef %>% magrittr::extract(,1),
                       se =~ t(.x)%*%vcv%*%.x %>% sqrt %>% magrittr::extract(,1)))) %>%
    rename_with(function(x){gsub("_inds", "", x)}, ends_with(c("se", "est"))) %>%
    select(-c(data, mod, vcv, coef, ends_with("inds"))) %>%
    pivot_longer(ends_with(c("est", "se")),
                 names_sep = "_",
                 names_to = c("par", ".value")) %>%
    mutate(pval = 2*pt(-abs(est/se), t_df)) %>%
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

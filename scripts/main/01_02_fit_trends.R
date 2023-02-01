source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

epa_regions <- read.csv(file.path(path_dropbox, "data", "us_climate_regions.csv"))
nonContig_stateFIPS <- c("02","60","66","15","72","78","69")

state_dict <- tigris::states(cb =  TRUE) %>% 
  filter(STATEFP %in% nonContig_stateFIPS == FALSE) %>% 
  left_join(epa_regions, by = c("STUSPS" = "state_code")) %>% 
  mutate(climate_regions = gsub(" \\(.*", "", climate_regions),
         climate_regions = stringr::str_to_title(climate_regions)) %>% 
  sf::st_drop_geometry() %>% 
  select(state_code = STATEFP, 
         state_abbr = STUSPS, 
         state_name = NAME,
         climate_regions = climate_regions)

station_year_avg <- readRDS(file.path(path_dropbox, "data", "epa_station_year_filled2022.rds")) %>% 
  rename(state_name = state) %>% 
  left_join(state_dict)

# break years from brandon ---- 
breaks <- rbind(read.csv(file.path(path_dropbox, "output", "year_subsample_estimates.csv")) %>%  
                  mutate(PM_measure = "mean"), 
                read.csv(file.path(path_dropbox, "output", "epm_year_subsample_estimates.csv")) %>% 
                  mutate(PM_measure = "extremePct")) %>%
  select(-X, -climate_regions) %>% 
  left_join(state_dict %>% 
              mutate(state_code = as.numeric(state_code))) %>% 
  select(-starts_with("state")) %>% 
  unique %>%
  rename_with(function(x) paste0("breakpoint-", x), 
              !c(climate_regions, PM_measure, starts_with("nrow"))) %>% 
  rename_with(function(x) gsub("nrow_", "nrow-", x), 
              starts_with("nrow")) %>% 
  pivot_longer(!c(climate_regions, PM_measure), 
               names_sep = "-", 
               names_to = c(".value", "spec"))

specs <- data.frame(spec = c("main", "break_1early", "break_1late", 
                             "allobs", "50obs5yr", "50obs15yr", "100obs10yr",
                             "drop2020", "drop2021", "drop2022",
                             "no_duplicated_break", "piecewise", 
                             "regional_breaks"),
                    split = c(0, -1, 1, 
                              0, 0, 0, 0,
                              0, 0, 0, 
                              0, 0, 
                              0), 
                    min_obs = c(50, 50, 50, 
                                1, 50, 50, 100, 
                                50, 50, 50, 
                                50, 50,
                                50), 
                    min_years = c(10, 10, 10, 
                                  1, 5, 15, 10, 
                                  10, 10, 10, 
                                  10, 10, 
                                  10), 
                    y_start = -Inf,
                    y_end = Inf,
                    y_drop = c(NA, NA, NA, 
                               NA, NA, NA, NA, 
                               2020, 2021, 2022, 
                               NA, NA, 
                               NA), 
                    fe_regex = c(rep("id|pm_type|period", 11), "id|pm_type", 
                                 "id|pm_type|period"),
                    duplicated_break = c(rep(TRUE, 10), 
                                         FALSE, FALSE,
                                         TRUE)) %>% 
  left_join(data.frame(pm_measure = c("mean", "extremePct"), 
                       main_break = c(2016, 2010)), 
            by = character()) %>% 
  mutate(split = main_break + split) %>% 
  select(-main_break)

all_fits <- purrr::pmap_dfr(
  specs, 
  function(spec, pm_measure, split, min_obs, min_years, 
           y_start, y_end, y_drop, fe_regex, 
           duplicated_break) {
    # if its not regional fit, just run fit_pooled onced
    if(spec != "regional_breaks"){
      spec_fit <- fit_pooled(station_year_avg, pm_measure, split, 
                             min_obs, min_years, 
                             y_start, y_end, y_drop, 
                             fe_regex, duplicated_break) %>% 
        mutate(breakpoint = split)
    } else{ # for regional breaks, run it for each region with that break point 
      spec_fit <- breaks %>% 
        filter(spec == "n50_y10" & PM_measure == pm_measure) %>% 
        transmute(region = climate_regions, 
                  region_breakpoint = round(breakpoint)) %>%
        purrr::pmap_dfr(function(region, region_breakpoint) {
          station_year_avg %>% 
            filter(climate_regions == region) %>% 
            fit_pooled(pm_measure, region_breakpoint, 
                       min_obs, min_years, 
                       y_start, y_end, y_drop, 
                       fe_regex, duplicated_break)  %>%
            mutate(breakpoint = region_breakpoint) %>% 
            filter(state_abbr != "US")
        })
    }
    spec_fit %>%
      mutate(spec = spec, 
             PM_measure = pm_measure) %>% 
      return
  }) %>% 
  pivot_wider(values_from = est:pval, names_from = par, names_glue = "{par}_{.value}") %>% 
  left_join(state_dict) %>%
  mutate(state_name = ifelse(state_abbr == "US", "Contiguous US", state_name))

# fits on the western states to assess sensitivity to breakpoints
west_fits <- purrr::map_dfr(
  2008:2016, 
  function(split) {
    station_year_avg %>% 
      filter(state_abbr %in% c("CA", "OR", "WA",  
                               # "MT", "UT", "CO", "WY",
                               "NV", "ID")) %>% 
      fit_pooled("mean", split, 50, 10, -Inf, Inf, NA, "id|pm_type|period", TRUE) %>%
      filter(state_abbr != "US") %>%
      mutate(split = split, 
             PM_measure = "mean") %>%
      return
  }) %>% 
  pivot_wider(values_from = est:pval, names_from = par, names_glue = "{par}_{.value}") %>% 
  left_join(state_dict)

# classify both and then save them 
west_fits %>% classify_trends() %>% 
  saveRDS(file.path(path_dropbox, "output", 
                    "west_breakpointFits_coefs_classifications.rds"))

all_fits %>% classify_trends() %>% 
  saveRDS(file.path(path_dropbox, "output", 
                    "state_allFits_coefs_classifications.rds"))

# ------------------------------------------------------------------------------
# Written by: Marissa Childs
# Calculates numbers in main text.
# ------------------------------------------------------------------------------
source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

# Set up ----
# sig level for distributions, two-sided
q_alpha = 0.05
main_spec = "n50_y15"

# Define color schemes 
smokeinfluence_colors <- c("#62205f", "#bb292c", "#ef8737", "#ffd353", "#dea868", "grey95")
extreme_colors <- c("grey90", "#5a5a83")
pmtrend_colors <- c("#175449", "#6e948c", "#c38f16", "#dcdcdb", "#5d6174")

# State dictionary ----
nonContig_stateFIPS <- c("02","60","66","15","72","78","69")
epa_regions <- read.csv(file.path(path_dropbox, "data", "us_climate_regions.csv"))

state_dict <- tigris::states(cb =  TRUE) %>% 
  filter(STATEFP %in% nonContig_stateFIPS == FALSE) %>% 
  left_join(epa_regions, by = c("STUSPS" = "state_code")) %>% 
  mutate(climate_regions = gsub(" \\(.*", "", climate_regions),
         climate_regions = stringr::str_to_title(climate_regions)) %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(state_code = STATEFP, 
                state_abbr = STUSPS, 
                state_name = NAME,
                climate_regions = climate_regions)

# Station-year data ----
station_year <- readRDS(file.path(path_dropbox, "data", "epa_station_year_2000_2022.rds")) %>% 
  rename(state_name = state) %>% 
  left_join(state_dict)

# State-year averages ----
state_year <- station_year %>% 
  filter(totalPM_n > 50) %>% 
  filter(n() >= 15, 
         .by = id) %>% 
  pivot_longer(contains("PM"),
               names_sep = "_",
               names_to = c("pm_type", ".value")) %>%
  {rbind(group_by(., state_abbr, state_name, climate_regions, year, pm_type) %>% 
           summarise(across(c(mean, extremePct), 
                            ~ mean(.x, na.rm = T)),
                     .groups = "drop"),
         group_by(., year, pm_type) %>% 
           summarise(across(c(mean, extremePct), 
                            ~ mean(.x, na.rm = T)),
                     .groups = "drop") %>%
           mutate(state_abbr = "US", 
                  state_name = "Contiguous US", 
                  climate_regions = NA))}

# Load bootstrap fits ----
mean_fits <- rbind(
  feather::read_feather(file.path(
    path_dropbox, "output", "conus", 
    "conus_mean_1000b_FErmperfect_alpha20_tol1e-08_corner95_fitstatdev_smartinit__maxit5_maxrestart2_COMBINED.feather")) %>% 
    mutate(boot = 1:n(), .by = c(sample_name, state)),
  feather::read_feather(file.path(
    path_dropbox, "output", "regional", 
    "regional_mean_1000b_FErmperfect_alpha15_tol1e-08_corner95_fitstatdev_smartinit__maxit5_maxrestart2_COMBINED.feather")) %>% 
    mutate(sample_name = paste0("regional_", sample_name)) %>% 
    mutate(boot = 1:n(), .by = c(sample_name, state)),
  feather::read_feather(file.path(
    path_dropbox, "output", "other", 
    "conus_discontinuous_mean_1000b_alpha15_tol1e-08_corner95_fitstatwll_smartinit__maxit5_maxrestart2_n50_y15.feather")) %>% 
    mutate(sample_name = "discontinuous", 
           it_algo = NA),
  feather::read_feather(file.path(
    path_dropbox, "output", "other", 
    "anom_compare_conus_mean_1000b_alpha20_corner95_fitstatdev_smartinit_maxit5_maxrestart2_pos_anom.feather")) %>% 
    mutate(boot = 1:n(), .by = c(sample_name, state))
) %>%
  add_sample_factor(main_spec = main_spec) %>%  
  # Filter out if either the total or nonsmoke regression didn't converge
  filter(if_all(starts_with("conv"), ~ .x) & if_all(matches("^b1|^b2"), ~!is.na(.x))) %>% 
  # Calculate coefficient differences
  mutate(b2nonsmokeMb2total = b2nonsmoke - b2total, 
         b2totalMb1total = b2total - b1total)

extreme_fits <- feather::read_feather(file.path(path_dropbox, "output", "conus", "conus_extremeN_1000b_FErmperfect_alpha20_tol1e-08_corner95_fitstatdev_smartinit__maxit5_maxrestart2_COMBINED.feather")) %>% 
  rbind(feather::read_feather(file.path(
    path_dropbox, "output", "other", 
    "anom_compare_conus_extremeN_1000b_alpha20_corner95_fitstatdev_smartinit_maxit5_maxrestart2_pos_anom.feather"))) %>% 
  add_sample_factor(main_spec = main_spec) %>% 
  # Recode to zero coef if no extreme days in the sample, 
  mutate(b2total = ifelse(mean_post_totalPM == 0, 0, b2total),
         b2nonsmoke = ifelse(mean_post_nonsmokePM == 0, 0, b2nonsmoke),
         # Set convergence = T since we only care about post-period slopes for 
         # the extreme fits
         conv_total = ifelse(mean_post_totalPM == 0, TRUE, conv_total),
         conv_nonsmoke = ifelse(mean_post_nonsmokePM == 0, TRUE, conv_nonsmoke)) %>%
  # Filter out if either the total or nonsmoke regression didn't converge
  filter(if_all(starts_with("conv"), ~ .x) & if_all(matches("^b2"), ~!is.na(.x))) %>% 
  # Calculate coefficient differences
  mutate(b2nonsmokeMb2total = b2nonsmoke - b2total)

# Classify fits ----
mean_classified <- mean_fits %>% 
  dplyr::select(state, sample_name, sample_factor, break_year, starts_with("b"), -boot) %>% 
  group_by(state, sample_name, sample_factor) %>%
  classify_boot(full = T)
extreme_classified <- extreme_fits %>% 
  dplyr::select(state, sample_name, sample_factor, break_year, starts_with("b2")) %>% 
  group_by(state, sample_name, sample_factor) %>%
  classify_boot(full = F)

mean_break <- mean_fits %>% 
  # Just so we only get 1 obs per bootstrap
  filter(state == "US") %>% 
  filter(sample_name == main_spec) %>%  
  pull(break_year) %>% median
extreme_break <- extreme_fits %>% 
  # Just so we only get 1 obs per bootstrap
  filter(state == "US") %>% 
  filter(sample_name == main_spec) %>% 
  pull(break_year) %>% median

# Calculate percent undone ----
pct_undone <- mean_fits %>%
  filter(state != "US" & sample_name == main_spec) %>% 
  mutate(boot = 1:n(), 
         .by = c(sample_name, state)) %>% 
  left_join(mean_classified %>% 
              filter(sample_name == main_spec) %>% 
              dplyr::select(state = state_name, 
                            smoke_influenced, smoke_group)) %>% 
  filter(smoke_influenced) %>% dplyr::select(-smoke_influenced) %>% 
  transmute(boot, smoke_group, state,
            changeEarlyTotal = b1total*(break_year - 2000),
            changeLateSmoke = (-b2nonsmokeMb2total)*(2022 - break_year),
            pctUndone = pmax(pmin(-changeLateSmoke/changeEarlyTotal, 1), 0)*100, 
            yearsUndone = changeLateSmoke/b1total)
pct_undone_groups <- pct_undone %>% 
  dplyr::select(-state) %>% 
  {rbind(., 
         mutate(., smoke_group = "all_smoke-influenced"), 
         filter(., grepl("reversal", smoke_group)) %>% 
           mutate(smoke_group = "smoke_reversals"))} %>% 
  # For each bootstrap, calculate avg (median) pct undone, 
  summarise(across(everything(), 
                   list(mid =~ median(.x), 
                        mean =~ mean(.x))), 
            .by = c(boot, smoke_group)) %>% 
  summarise(across(!boot, 
                   list(lower =~ quantile(.x, q_alpha/2), 
                        mid =~ quantile(.x, 0.5),
                        # mean =~ quantile(.x, 0.5),
                        upper =~ quantile(.x, 1 - q_alpha/2))), 
            .by = smoke_group) %>% 
  pivot_longer(!smoke_group) %>% 
  separate_wider_delim(cols = name, 
                       names = c("name", "sample_stat", "boot_stat"), 
                       delim = "_") %>% 
  mutate(panel = ifelse(grepl("^change", name), 
                        "change in PM2.5",
                        "pct of progress undone by smoke"))

# Load station-day data ----
station_day <- readRDS(file.path(path_dropbox, "data", "epa_station_day_totalPM_smokePM_panel_20000101-20221231.rds")) %>% 
  filter(!is.na(pm25))

# Percent of exceedances from wildfire smoke (early vs late on breakpoint) and 
# post 2020 (last 3 years) ----
extremePct_cdf <- station_year %>% 
  filter(totalPM_n > 50) %>% 
  filter(n() >= 15, 
         .by = id) %>% 
  mutate(pct_exceed = ifelse(totalPM_extremePct > 0, 
                             (totalPM_extremePct - nonsmokePM_extremePct)/totalPM_extremePct, 
                             0)) %>% 
  {purrr::pmap_dfr(data.frame(start = c(2006, 2011, 2020, 2006:2022),
                              end = c(2010, 2022, 2022, 2006:2022)),
                   function(start, end){
                     filter(., year >= start & year <= end) %>%
                       group_by(year, state_abbr) %>% 
                       summarise(across(ends_with("exceed"), 
                                        mean),
                                 .groups = "drop") %>% 
                       group_by(state_abbr) %>% 
                       summarise(across(ends_with("exceed"), 
                                        mean),
                                 .groups = "drop") %>% 
                       arrange(desc(pct_exceed)) %>% 
                       mutate(n_gte = 1:n()) %>% 
                       group_by(pct_exceed) %>%
                       summarise(n_gte = max(n_gte), 
                                 states = list(state_abbr)) %>%
                       # Add row to get down to zero 
                       rbind(., 
                             filter(., pct_exceed == max(pct_exceed)) %>% 
                               mutate(n_gte = 0)) %>%
                       mutate(start = start, 
                              end = end) %>% 
                       return
                   })} %>% 
  
  mutate(period = ifelse(start == end, 
                         start, 
                         paste0(start, "-", end)),
         type = ifelse(start == end, "annual", "range"),
         group = ifelse(start == end, period, "range")) 

# Load vpd data ----
historical_vpd <- read_rds(file.path(path_dropbox, 'data', 
                                     'gridmet_vpd_1979_2022_1deg_forest_cover.rds')) %>% 
  filter(west == 1)
vpd_obs_month <- historical_vpd %>% 
  filter(month > 4, month < 10) %>% 
  group_by(year, month) %>% 
  summarise(vpd = weighted.mean(vpd, forest, na.rm = T)) %>% 
  ungroup()
vpd_obs_year <- vpd_obs_month %>% 
  group_by(year) %>% 
  summarise(vpd = mean(vpd, na.rm = T)) %>% 
  ungroup() 

# Load smokePM ----
dt <- read_rds(file.path(path_dropbox, "data", "epa_station_year_2000_2022.rds"))

# Define Western states ----
ws <- c("California", "Oregon", "Washington", "Nevada", "Arizona", "New Mexico", 
        "Utah", "Idaho", "Montana", "Wyoming", "Colorado")
dty <- dt %>% 
  mutate(smokePM = totalPM_mean-nonsmokePM_mean) %>% 
  filter(state %in% ws) %>% 
  group_by(year) %>% 
  summarise(smokePM = mean(smokePM,na.rm = T))
tp = full_join(dty, vpd_obs_year)

# Load gcm data ----
gcm <- read_rds(file.path(path_dropbox, 'data', 'CMIP6_west_vpd_debias_annual_2015_2100.rds'))
gcms <- gcm %>% 
  filter(year >= 2040 & year <= 2060) %>%
  group_by(model, scenario) %>% 
  summarise(vpd = mean(vpd_debias))

# Compute change in smoke at ensemble median VPD ----
ss <- c("ssp126", "ssp245", "ssp370")
meds <- gcms %>% 
  filter(scenario %in% ss) %>% 
  group_by(scenario) %>% 
  summarise(vpd = median(vpd))
meds <- rbind(meds, 
              data.frame(scenario = "base", 
                         vpd = mean(tp$vpd[tp$year >= 2018])))

summary(mod <- lm(log(smokePM) ~ vpd, tp))
s = summary(mod)$sigma
a0 = sum(exp(mod$residuals))/length(mod$residuals) # smearing estimate, eq 6.43 in wooldridge
meds <- meds %>% 
  mutate(pred = exp(predict(mod, meds))*exp(s^2/2), 
         predsmear = a0*exp(predict(mod, meds)))

data.frame(scen = meds$scenario[1:3], predsmear = meds$predsmear[1:3] - meds$predsmear[4])

gcmp <- data.frame(gcms) %>% 
  mutate(pred = a0*exp(predict(mod, data.frame(vpd = gcms$vpd)))) %>% 
  filter(scenario != "ssp585") %>% 
  group_by(scenario) %>% 
  summarise(meanpred = mean(pred)) %>% 
  left_join(., meds)

# ------------------------------------------------------------------------------
# main text numbers ----
out_txt_file <- file.path(path_github, "tables", "raw", "manuscript_numbers.txt")


{# pct of states, smoke-influence + stagnating/reversing 
  write("abstract, % of states where smoke influence is slowing/reversing progress ----------", 
        out_txt_file)
  mean_classified %>% 
    filter(sample_name == main_spec & state_abbr != "US") %>% 
    # mutate(group = ifelse(grepl("no", smoke_group_factor), "no influence/no trend impact", "smoke-influenced, stagnating or reversing")) %>%
    summarise(pct = n()/48, 
              .by = smoke_influenced) %>% 
    write.table(out_txt_file,
                append = TRUE, row.names = FALSE) 
  
  # abstract, pct gains undone 
  write("\n\n\nabstract, % of gains undone ----------", 
        out_txt_file, 
        append = T)
  pct_undone_groups %>% 
    filter(name == "pctUndone" & sample_stat == "mean" & 
             boot_stat == "mid" & smoke_group == "all_smoke-influenced") %>% 
    pull(value) %>% 
    as.character() %>% 
    write(out_txt_file, 
          append = T)
  
  # abstract, years undone
  write("\n\n\nabstract, mean years undone by smoke over all smoke-influenced states ----------", 
        out_txt_file, 
        append = T)
  pct_undone_groups %>% 
    filter(name == "yearsUndone" & sample_stat == "mean" & 
             boot_stat == "mid" & smoke_group == "all_smoke-influenced") %>% 
    pull(value) %>% 
    as.character() %>% 
    write(out_txt_file, 
          append = T)
  
  # intro, break years
  write("\n\n\nintro, median break year for annual average PM, main specification ----------", 
        out_txt_file, 
        append = T)
  mean_fits %>% 
    filter(state == "US" & sample_name == main_spec) %>% 
    pull(break_year) %>% 
    median %>% 
    format(digits = 12) %>% 
    write(out_txt_file, 
          append = T)
  
  write("\n\n\nintro, median break year for extreme PM days, main specification ----------", 
        out_txt_file, 
        append = T)
  extreme_fits %>% 
    filter(state == "US" & sample_name == main_spec) %>% 
    pull(break_year) %>% 
    median %>% 
    format(digits = 12) %>% 
    write(out_txt_file, 
          append = T)
  
  # intro, eastern states extreme days 
  # number of extreme days in 2000 - 2002 in eastern states 
  write("\n\n\nintro, extreme days in eastern states 2000-2002 ----------", 
        out_txt_file, 
        append = T)
  station_year %>% 
    filter(year %in% 2000:2002) %>% 
    filter(climate_regions %in% c("Northeast", "Ohio Valley")) %>%
    summarise(avg_extremePct = mean(totalPM_extremePct)) %>% 
    multiply_by(365) %>% 
    magrittr::extract(1, 1) %>% 
    write(out_txt_file, 
          append = T)
  
  # results, counts of total PM trends 
  write("\n\n\nresults, counts of states with PM trend changes ----------", 
        out_txt_file, 
        append = T)
  mean_classified %>% 
    filter(sample_name == main_spec & state_abbr != "US") %>% 
    summarise(n = n(), 
              .by = change_group) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  # results, pct smoke influenced
  write("\n\n\nresults, pct and number of smoke-influenced states ----------", 
        out_txt_file, 
        append = T)
  mean_classified %>% 
    filter(sample_name == main_spec & state_abbr != "US") %>% 
    summarise(pct = n()/48, 
              n = n(),
              .by = smoke_influenced) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  # results, PM added from smoke 
  write("\n\n\nresults, ug/m3 PM added from smoke in smoke-influenced states ----------", 
        out_txt_file, 
        append = T)
  pct_undone_groups %>% 
    filter(name == "changeLateSmoke" & smoke_group == "all_smoke-influenced") %>% 
    mutate(value = round(value, 2)) %>% 
    pivot_wider(names_from = boot_stat, values_from = value) %>% 
    # select(-smoke_group, -name, -panel) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  # results, pct undone
  write("\n\n\nresults, pct progress undone in all smoke-influenced states ----------", 
        out_txt_file, 
        append = T)
  pct_undone_groups %>% 
    filter(name == "pctUndone" & smoke_group == "all_smoke-influenced") %>% 
    # mutate(value = round(value, 0)) %>% 
    pivot_wider(names_from = boot_stat, values_from = value) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  # results, PM added and pct undone among smoke-influenced groups
  write("\n\n\nresults, recap state counts by smoke-influenced group ----------", 
        out_txt_file, 
        append = T)
  mean_classified %>% 
    filter(sample_name == main_spec & state_abbr != "US") %>% 
    summarise(n = n(),
              .by = smoke_group) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  write("\n\n\nresults, PM added and pct progress undone by smoke group ----------", 
        out_txt_file, 
        append = T)
  pct_undone_groups %>% 
    filter(smoke_group %in% c("smoke-influenced stagnation", 
                              "smoke_reversals") & 
             sample_stat == "mid" &
             name %in% c("changeLateSmoke", 
                         "pctUndone", 
                         "yearsUndone")) %>% 
    mutate(value = round(value, case_when(name == "pctUndone" ~ 0,
                                          name == "changeLateSmoke" ~ 2,
                                          name == "yearsUndone" ~ 1))) %>% 
    pivot_wider(names_from = boot_stat, values_from = value) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  # results, smoke-influenced sensitivity 
  write("\n\n\nresults, smoke-influenced group sensitivity ----------", 
        out_txt_file, 
        append = T)
  mean_classified %>% 
    filter(sample_name %in% c("discontinuous", "drop2021") & 
             state_abbr != "US") %>% 
    summarise(n = n(), 
              .by = c(sample_name, smoke_influenced)) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  # results, extremes smoke-influenced 
  write("\n\n\nresults, extreme days smoke influenced counts ----------", 
        out_txt_file, 
        append = T)
  extreme_classified %>% 
    filter(sample_name == main_spec & state_abbr != "US") %>% 
    summarise(n = n(), 
              .by = c(smoke_influenced)) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  # results, extremes smoke-influenced sensitivity 
  write("\n\n\nresults, extreme days sensitivity ----------", 
        out_txt_file, 
        append = T)
  extreme_classified %>% 
    summarise(n = n(), 
              .by = c(sample_factor, smoke_influenced)) %>% 
    filter(smoke_influenced == TRUE) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  # results, extreme pct numbers 
  write("\n\n\nresults, pct of states with >25% exceedences from smoke by time period  ----------", 
        out_txt_file, 
        append = T)
  extremePct_cdf %>% 
    filter(pct_exceed >= 0.25) %>% 
    group_by(period) %>% 
    summarize(n = max(n_gte), 
              states = list(unique(states))) %>% 
    filter(grepl("-", period)) %>%
    unnest_wider(states, names_sep = "") %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  write("\n\n\nresults, states with >75% exceedences from smoke in 2020-2022 ----------", 
        out_txt_file, 
        append = T)
  extremePct_cdf %>% 
    filter(pct_exceed >= 0.75) %>% 
    filter(period == "2020-2022") %>%
    summarize(n = max(n_gte), 
              states = list(unique(states))) %>% 
    unnest_wider(states, names_sep = "") %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE)
  
  # results, trend in western states  
  write("\n\n\nresults, total PM trend in western states ----------", 
        out_txt_file, 
        append = T)
  feols(totalPM_mean ~ years_pre + years_post | id, 
        data = station_year %>% 
          filter(totalPM_n > 50) %>% 
          filter(n() >= 15, 
                 .by = id) %>% 
          filter(state_abbr %in% c("CA", "OR", "WA", 
                                   "ID", "NV", "UT", 
                                   "AZ", "MT", "WY",
                                   "CO", "NM")) %>% 
          mutate(period = (year > mean_break),
                 years_post = (year - mean_break)*(year > mean_break), 
                 years_pre = (year - mean_break)*(year <= mean_break))) %>% 
    coef() %>% 
    magrittr::extract("years_post") %>% 
    write(out_txt_file, 
          append = T)
  
  
  # methods numbers 
  write("\n\n\nmethods, number of stations ----------", 
        out_txt_file, 
        append = T)
  station_year$id %>% unique %>% length %>% 
    as.character() %>% 
    write(out_txt_file, 
          append = T)
  
  write("\n\n\nmethods, pct of total days with negative smoke anomalies ----------", 
        out_txt_file, 
        append = T)
  station_day %>% 
    filter(!is.na(pm25_anom)) %>% 
    mutate(neg_smokePM = (smoke_day*pm25_anom) < 0) %>% 
    summarise(pct_total = (sum(neg_smokePM)/n())*100) %>% 
    magrittr::extract(1,1) %>% 
    round(2) %>% 
    paste0("%") %>% 
    write(out_txt_file, 
          append = T)
  
  write("\n\n\nmethods, avg increase in PM2.5 when smoke plume is overhead ----------", 
        out_txt_file, 
        append = T)
  feols(pm25 ~ smoke_day | id + state^month + state^year, 
        data = station_day %>% 
          filter(!is.na(smoke_day)) %>% 
          mutate(year = lubridate::year(date),
                 month = lubridate::month(date))) %>% 
    coeftable %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = TRUE)
  
  write("\n\n\nmethods, change in smoke-influenced classifications in other specifications ----------", 
        out_txt_file, 
        append = T)
  rbind(extreme_classified %>% 
          filter(state_abbr != "US") %>% 
          transmute(type = "extreme", 
                    smoke_influenced_text, 
                    sample_factor),
        mean_classified %>% 
          filter(state_abbr != "US") %>% 
          transmute(type = "mean", 
                    smoke_influenced_text, 
                    sample_factor)) %>% 
    summarise(n = n(), 
              .by = c(smoke_influenced_text, sample_factor, type)) %>% 
    pivot_wider(names_from = smoke_influenced_text, values_from = n, 
                values_fill = 0) %>% 
    arrange(desc(type), `smoke-influenced`) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = FALSE) 
  
  write("\n\nmethods, beta in regression of log smoke pm on vpd ----------", 
        out_txt_file, 
        append = T)
  # without time trend
  lm(log(smokePM) ~ vpd, tp) %>% 
    summary() %>% 
    coeftable() %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = TRUE)
  # with time trend
  lm(log(smokePM) ~ vpd + year, tp) %>% 
    summary() %>% 
    coeftable() %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = TRUE)
  
  write("\n\nmethods, predictions for future annual avg smoke pm under SSPs ----------", 
        out_txt_file, 
        append = T)
  gcmp %>% 
    select(scenario, predsmear) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = TRUE)
  
  write("\n\nmethods, predicted annual avg smoke pm from average VPD 2018-2022, and difference ----------", 
        out_txt_file, 
        append = T)
  # predicted annual avg smoke pm from average VPD 2018-2022
  meds %>% 
    filter(scenario == "base") %>% 
    pull(predsmear) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = TRUE)
  # difference
  ((gcmp %>% 
      filter(scenario == "ssp370") %>% 
      pull(predsmear)) - 
      (meds %>% 
         filter(scenario == "base") %>% 
         pull(predsmear))
  ) %>% 
    write.table(out_txt_file, 
                append = T, 
                row.names = TRUE)
}

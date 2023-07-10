# ------------------------------------------------------------------------------
# Written by: Marissa Childs
# Plots Figure ED1.
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

# Get station locations----
station_ll <- read_sf(file.path(path_dropbox, "data", "epa_station_locations"))

states <- tigris::states(cb = TRUE, year = 2020) %>% 
  filter(STATEFP %in% nonContig_stateFIPS == FALSE)

# ------------------------------------------------------------------------------
# Figure ED1: smoke PM definition and station locations ----
# Panel with daily station data 
station_day <- readRDS(file.path(path_dropbox, "data", 
                                 "epa_station_day_totalPM_smokePM_panel_20000101-20221231.rds")) %>% 
  filter(!is.na(pm25))
station_day %>% 
  filter(id == "060811001") %>%
  filter(date >= as.Date("2020-07-01") & date <= as.Date("2020-10-02")) %>%
  mutate(totalPM = pm25,
         nonsmokePM = totalPM - smoke_day*pm25_anom) %>% 
  dplyr::select(date, totalPM, nonsmokePM, smoke_day) %>% 
  mutate(nonsmoke_max = nonsmokePM, 
         nonsmoke_min = 0, 
         total_max = totalPM, 
         total_min = nonsmokePM) %>% 
  pivot_longer(ends_with(c('max', 'min')), 
               names_sep = "_", 
               names_to = c("pm_type", ".value")) %>%  
  mutate(pm_type = gsub("non", "non-", pm_type)) %>%
  arrange(pm_type) %>% 
  {ggplot(data = ., 
          aes(x = date)) + 
      geom_line(aes(y = max, color = pm_type, group = pm_type)) +
      geom_ribbon(aes(ymax = max, ymin = min, fill = pm_type), 
                  color = NA, alpha = 0.7) + 
      geom_point(aes(x = date, y = -3, 
                     alpha = I(ifelse(smoke_day, 1, 0))), 
                 color = "gray") + 
      scale_color_manual(values = c("blue", "black"),
                         aesthetics = c("fill", "color")) + 
      theme_classic() + 
      theme(legend.position = "none") + 
      xlab("") + 
      ylab(expression(paste("daily ", PM[2.5]," (",mu, "g/", m^3, ")")))} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureED01a_station-day_example.pdf"), 
         ., width = 3.5, height = 2.5)

# Panel with station-year non-smoke and smoke averages
station_year %>% 
  filter(id == "060811001") %>%
  pivot_longer(ends_with("mean")) %>% 
  mutate(name = gsub("PM_mean", "", name)) %>% 
  arrange(rev(name)) %>% 
  {ggplot(data = ., 
          aes(x = year, y = value, color = name, 
              linetype = name)) + 
      geom_line(lwd = 0.8) + 
      scale_color_manual(values = c("blue", "black"),
                         aesthetics = c("fill", "color")) + 
      scale_linetype_manual(values = c("32", "solid")) +
      theme_classic() + 
      theme(legend.position = "none") + 
      xlab("") + 
      ylab(expression(paste("annual ", PM[2.5]," (",mu, "g/", m^3, ")")))} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureED01b_station-year_example.pdf"), 
         ., width = 3.5, height = 2.5)

# Panel with station locations by number of years with > 50 obs
station_year %>% 
  filter(totalPM_n > 50) %>% 
  group_by(id) %>% 
  summarise(n_year = n()) %>% 
  mutate(n_year = case_when(n_year < 5 ~ "< 5 years", 
                            n_year >= 5 & n_year < 10 ~ "5 - 9 years", 
                            n_year >= 10 & n_year < 15 ~ "10 - 14 years", 
                            n_year >= 15 ~ "15+ years"), 
         n_year = factor(n_year, 
                         levels = c("15+ years", 
                                    "10 - 14 years", 
                                    "5 - 9 years", 
                                    "< 5 years"), 
                         ordered = TRUE)) %>%
  left_join(station_ll %>% dplyr::select(id = stn_id)) %>% 
  st_as_sf() %>% 
  st_transform(st_crs(states)) %>% 
  arrange(desc(n_year)) %>% 
  {ggplot(data = .) + 
      geom_sf(data = states, fill = NA) +
      geom_sf(aes(color = n_year), size = 0.8, alpha = 0.8) +
      scale_color_manual(name = "years with\nat least\n50 observations",
                         values = c("#591c19", "#b64f32", "#f7c267", "#8b8b99")) +
      theme_void() + 
      theme(legend.position = "right")} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureED01c_station_locations.pdf"), 
         ., width = 7, height = 4)

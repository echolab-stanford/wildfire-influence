source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

# set up state grid
my_grid <- us_state_grid3
# move ME up, VT and NH one right, FL one left
my_grid %<>% 
  mutate(row = ifelse(code == "ME", row - 1, row), 
         col = ifelse(code %in% c("VT", "NH"), col + 1, col), 
         col = ifelse(code == "FL", col - 1, col)) %>% 
  # remove AK and HI
  filter(code %in% c("AK", "HI") == FALSE) %>% 
  # add US to top left
  rbind(data.frame(row = 1, col = 1, code = "US", name = "Contiguous US"))

# state dictionary ----
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

# station-year data ----
station_year <- readRDS(file.path(path_dropbox, "data", 
                                  "epa_station_year_filled2022.rds")) %>% 
  rename(state_name = state) %>% 
  left_join(state_dict)

# state-year averages ----
state_year <- station_year %>% 
  filter(totalPM_n > 50) %>% 
  group_by(id) %>% 
  mutate(n_year = n()) %>% 
  ungroup %>% 
  filter(n_year >= 15) %>% 
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

# read in fits ----
all_fits <- readRDS(file.path(path_dropbox, "output", 
                              "state_allFits_coefs_classifications.rds")) %>% 
  mutate(smoke_group = gsub("no smoke-influence detected", 
                            "no smoke influence detected", 
                            smoke_group)) %>% 
  mutate(smoke_group = gsub(", ", ",\n", smoke_group),
         smoke_group_factor = factor(smoke_group, 
                                     levels = c("smoke-caused reversal", 
                                                "smoke-influenced reversal", 
                                                "smoke-influenced stagnation", 
                                                "smoke-influenced,\nno early decline", 
                                                "smoke-influenced,\nno signficant trend change", 
                                                "no smoke influence detected"), 
                                     ordered = TRUE), 
         change_group = gsub("sig", "sig.", change_group),
         change_group_factor = factor(change_group, 
                                      levels = c("reversal", 
                                                 "stagnation",
                                                 "acceleration", 
                                                 "non-sig. change", 
                                                 "no sig. early decline"), 
                                      ordered = TRUE),
         smoke_influenced_text = ifelse(smoke_influenced, 
                                        "smoke-influenced", 
                                        "no smoke influence detected"), 
         spec_factor = ifelse(spec %in% c("break_1early", "break_1late") & 
                                PM_measure == "extremePct", 
                              paste0(spec, "_", PM_measure), 
                              spec), 
         spec_factor = recode_factor(spec_factor, 
                                     "main" = "Main", 
                                     "allobs" = "All obs", 
                                     "50obs5yr" = "5yrs, 50obs",
                                     "50obs10yr" = "10yrs, 50obs",
                                     "50obs15yr" = "15yrs, 50obs", 
                                     "100obs5yr" = "5yrs, 100obs", 
                                     "100obs10yr" = "10yrs, 100obs", 
                                     "break_1early" = "break in 2015", 
                                     "break_1late" = "break in 2017", 
                                     "break_1early_extremePct" = "break in 2009", 
                                     "break_1late_extremePct" = "break in 2011", 
                                     "drop2020" = "drop 2020", 
                                     "drop2021" = "drop 2021", 
                                     "drop2022" = "drop 2022", 
                                     "no_duplicated_break" = "no dup. break yr", 
                                     "piecewise" = "piecewise", 
                                     "regional_breaks" = "regional breaks",
                                     .ordered = TRUE))

# fig S1: smoke PM definition and station locations ----
# one panel with daily station data 
station_day <- readRDS(file.path(
  path_dropbox, "data", 
  "epa_station_day_totalPM_smokePM_panel_20000101-20221021_era5.rds"
)) %>% 
  filter(!is.na(pm25)) 
station_day %>% 
  filter(id == "060811001") %>%
  filter(date >= as.Date("2020-07-01") & date <= as.Date("2020-10-02")) %>%
  mutate(nonsmokePM = pm25 - smokePM, 
         totalPM = pm25) %>% 
  select(date, totalPM, nonsmokePM, smoke_day) %>% 
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
  ggsave(file.path(path_github, "figures", "raw", "figureS01a_station-day_example.pdf"), 
         ., width = 3.5, height = 2.5)

# second with station-year non-smoke and smoke averages
station_year %>% 
  # filter(id == "060010009") %>%
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
  ggsave(file.path(path_github, "figures", "raw", "figureS01b_station-year_example.pdf"), 
         ., width = 3.5, height = 2.5)

# lower panel with station locations by number of years with > 50 obs
station_ll <- read_sf(file.path(path_dropbox, "data", "epa_station_locations"))

states <- tigris::states(cb = TRUE, year = 2020) %>% 
  filter(STATEFP %in% nonContig_stateFIPS == FALSE)

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
  left_join(station_ll %>% select(id = stn_id)) %>% 
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
  ggsave(file.path(path_github, "figures", "raw", "figureS01c_station_locations.pdf"), 
         ., width = 7, height = 4)

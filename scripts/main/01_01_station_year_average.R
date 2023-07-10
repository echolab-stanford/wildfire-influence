# ------------------------------------------------------------------------------
# Written by: Marissa Childs
# Aggregates smoke PM2.5 to the annual level.
# ------------------------------------------------------------------------------
source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

# ------------------------------------------------------------------------------
# Load daily station data
station_day <- readRDS(file.path(
  path_dropbox, "data", 
  "epa_station_day_totalPM_smokePM_panel_20000101-20221231.rds"
))

# Calculate non-smoke assuming only positive pm anomalies are smoke PM 
station_day_pos <- station_day %<>%
  filter(!is.na(pm25)) %>%
  mutate(year = lubridate::year(date),
         totalPM = pm25,
         nonsmokePM = totalPM - smoke_day*pmax(pm25_anom, 0))

# Calculate non-smoke assuming ALL pm anomalies are smoke PM (even negative anomalies)
station_day %<>%
  filter(!is.na(pm25)) %>%
  mutate(year = lubridate::year(date),
         totalPM = pm25,
         nonsmokePM = totalPM - smoke_day*pm25_anom)

# Calculate station-year averages for both positive and all anomalies 
station_year_avg <- station_day %>% 
  select(state_code, state, id, year, totalPM, nonsmokePM) %>% 
  group_by(state_code, state, id, year) %>%
  # Calculate the average annual PM2.5 and the pct of days that are extreme (> 35ug)
  summarise(across(everything(), 
                   list(mean =~ mean(.x, na.rm = T), 
                        extremePct =~ sum(.x >= 35, na.rm = T)/length(na.omit(.x)),
                        n =~ length(na.omit(.x)))), 
            .groups = "drop")

station_year_avg_pos <- station_day_pos %>% 
  select(state_code, state, id, year, totalPM, nonsmokePM) %>% 
  group_by(state_code, state, id, year) %>%
  summarise(across(everything(), 
                   list(mean =~ mean(.x, na.rm = T), 
                        extremePct =~ sum(.x >= 35, na.rm = T)/length(na.omit(.x)),
                        n =~ length(na.omit(.x)))), 
            .groups = "drop")

# Save data 
saveRDS(station_year_avg_pos,
        file.path(path_dropbox, "data", "epa_station_year_2000_2022_posAnom.rds"))

saveRDS(station_year_avg,
        file.path(path_dropbox, "data", "epa_station_year_2000_2022.rds"))

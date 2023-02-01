source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

station_day <- readRDS(file.path(
  path_dropbox, "data", 
  "epa_station_day_totalPM_smokePM_panel_20000101-20221021_era5.rds"
)) %>% 
  filter(!is.na(pm25)) %>%
  mutate(year = lubridate::year(date),
         totalPM = pm25,
         nonsmokePM = totalPM - smokePM) 

# calculate station-year average
station_year_avg <- station_day %>% 
  select(state_code, state, id, year, totalPM, nonsmokePM) %>% # , temperature:wind_v_850hpa
  group_by(state_code, state, id, year) %>%
  summarise(across(everything(), 
                   list(mean =~ mean(.x, na.rm = T), 
                        extremePct =~ sum(.x >= 35, na.rm = T)/length(na.omit(.x)),
                        n =~ length(na.omit(.x)))), 
            .groups = "drop")

# for the last year, find the date range of the missings obs, for each station 
# for the last 3 years, average PM2.5 and smokePM2.5 over those days, then 
# replace the values
test <- station_day %>% 
  group_by(year) %>% 
  summarise(start = min(date), 
            end = max(date)) %>% 
  mutate(doy = lubridate::yday(end)) %>% 
  filter(doy != 365 + lubridate::leap_year(year))
print(test)

# want to average oct 22 - dec 31 for the last 3 years, accounting for leap year
station_fill <- station_day %>%
  mutate(doy_noleap = yday(date), 
         leap_correct = leap_year(date)*(doy_noleap >= yday(as.Date("2000-03-01"))), 
         doy_noleap = doy_noleap + ifelse(leap_correct, -1, 0)) %>% 
  # this isn't currently robust to the year with missing data being a leap year
  filter(doy_noleap > test$doy[1] & year %in% 2019:2021) %>% 
  select(state_code, id, year, totalPM, nonsmokePM) %>% # , temperature:wind_v_850hpa
  group_by(state_code, id, year) %>%
  summarise(across(everything(), 
                   list(mean =~ mean(.x, na.rm = T),
                        extremePct =~ sum(.x >= 35, na.rm = T)/length(na.omit(.x)),
                        n =~ length(na.omit(.x)))), 
            .groups = "drop")  %>% 
  group_by(state_code, id) %>% 
  summarise(across(contains("PM"), 
                   ~mean(.x, na.rm = T)), 
            .groups = "drop")

# should update this so it isn't based on manually calculating the numbers, 
# but i think eventually we'll have full 2022 data, so less of an issue
station_year_avg %<>% 
  left_join(station_fill %>% 
              rename_with(~paste0(.x, "_fill"), 
                          contains("PM")) %>% 
              mutate(year = 2022)) %>% 
  mutate(totalPM_mean = ifelse(year == 2022, 294/365*totalPM_mean + 71/365*totalPM_mean_fill, totalPM_mean), 
         nonsmokePM_mean = ifelse(year == 2022, 294/365*nonsmokePM_mean + 71/365*nonsmokePM_mean_fill, nonsmokePM_mean), 
         totalPM_n = ifelse(year == 2022, totalPM_n + totalPM_n_fill, totalPM_n), 
         nonsmokePM_n = ifelse(year == 2022, nonsmokePM_n + nonsmokePM_n_fill, nonsmokePM_n),
         totalPM_extremePct = ifelse(year == 2022, 294/365*totalPM_extremePct + 71/365*totalPM_extremePct_fill, totalPM_extremePct), 
         nonsmokePM_extremePct = ifelse(year == 2022, 294/365*nonsmokePM_extremePct + 71/365*nonsmokePM_extremePct_fill, nonsmokePM_extremePct)) %>% 
  select(-ends_with("fill")) #%>% 
# rename_with(function(x) {gsub("_mean", "", x)}, ends_with("mean"))

# save
saveRDS(station_year_avg,
        file.path(path_dropbox, "data", "epa_station_year_filled2022.rds"))

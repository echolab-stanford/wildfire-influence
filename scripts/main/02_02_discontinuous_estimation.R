# ------------------------------------------------------------------------------
# Written by: Brandon de la Cuesta
# Conducts discontinuous estimation.
# ------------------------------------------------------------------------------
# Set conus = T for US-wide breakpoint and F for regional
conus = T

# Set pm_type = "mean" for mean PM and "extremeN" (or anything not "mean") for extreme PM
# Only "mean" for discontinuous results
pm_type = "mean"

source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

# ------------------------------------------------------------------------------
# Load in data
df_epa <- read_csv(file.path(path_dropbox, "data", "us_climate_regions.csv")) %>% 
  mutate(
    state_code = fips(state_code), 
    climate_regions = case_when(
      climate_regions == "NORTHEN ROCKIES AND PLAINS (WEST NORTH CENTRAL)" ~ "WEST NORTH CENTRAL",
      climate_regions == "UPPER MIDWEST (EAST NORTH CENTRAL)" ~ "EAST NORTH CENTRAL",
      climate_regions == "OHIO VALLEY (CENTRAL)" ~ "CENTRAL",
      .default = climate_regions
    )
  ) %>% 
  mutate(climate_regions = tolower(gsub(" ", "_", climate_regions)))

df_year <- read_rds(file.path(path_dropbox, "data", 
                              "epa_station_year_2000_2022.rds")) %>%
  mutate(id_year = paste0(id, "-", year), 
         state_year = paste0(state, "-", year))
df_year <- left_join(df_year, df_epa, by = "state_code") %>%
  group_by(id) %>% 
  mutate(n_years = n_distinct(year), 
         totalPM_extremeN = totalPM_extremePct*totalPM_n,
         nonsmokePM_extremeN = nonsmokePM_extremePct*nonsmokePM_n)

# This dataset bottom codes smoke anomalies to 0
df_anom <- read_rds(file.path(path_dropbox, "data", 
                              "epa_station_year_2000_2022_posAnom.rds")) %>%
  mutate(id_year = paste0(id, "-", year), 
         state_year = paste0(state, "-", year))
df_anom <- left_join(df_anom, df_epa, by = "state_code") %>%
  group_by(id) %>% 
  mutate(n_years = n_distinct(year), 
         totalPM_extremeN = totalPM_extremePct*totalPM_n,
         nonsmokePM_extremeN = nonsmokePM_extremePct*nonsmokePM_n)

# Use cores_prop% of all cores
cores_prop <- 0.70
n_cores <- round(cores_prop*availableCores())
parallel <- T

if (pm_type != "mean") {stop("Only mean PM is supported for discontinuous fits.")}
print(rep(paste0("CURRENTLY DOING: DISCONTINUOUS ", ifelse(conus, "CONUS", "REGIONAL"), "_", pm_type), 5))

# Condition check to set other globals/global fixest options
if (parallel) {
  # If parallelizing, prevent fixest from using more than 1 thread at a time
  setFixest_nthreads(1)
  # Reset connections for each subsample to prevent dead workers
  closeAllConnections()
  registerDoFuture()
  plan(multisession, workers = n_cores)
}
# If not parallelizing, use cores_prop% of cores
if (!parallel) {setFixest_nthreads(n_cores)}

# Overwrite globals where necessary
# Use "wll" for both mean and extreme PM since it partials out the FE
break_fitstat <- "wll"
# For discontinuous models, should only be femlm
use_fun <- femlm
break_alpha <- 0.15
# Only main sample will be used here
sample_names <- "n50_y15"
comb_name <- file.path(
  path_dropbox, "output", 
  paste0(ifelse(
    conus, "other/conus_", "other/regional_"), 
    "discontinuous_", pm_type, "_", n_boot, "b", 
    "_alpha", break_alpha*100, "_tol", tol, "_corner", corner_push*100,  
    "_fitstat", break_fitstat, ifelse(smart_init, "_smartinit_", "_randominit_"), 
    "_maxit", max_it, "_maxrestart", max_restart, "_", sample_names, ".feather"))

# ------------------------------------------------------------------------------
#### Subset and declare formula ####
df_tmp <-  df_year %>% 
  filter(totalPM_n > 50) %>%
  group_by(id) %>% 
  mutate(n_years = n_distinct(year)) %>% 
  ungroup() %>%
  filter(n_years >= 15)

# Create station-level dataset, required for bootstrap function
df_station <- df_tmp %>% 
  group_by(id, state) %>% 
  summarise()

# State-level station counts
strat_n <- tapply(df_station$id, df_station$state, 
                  function(x) {length(unique(x))})

# ------------------------------------------------------------------------------
#### Estimation #####
df_full <- foreach(i = 1:n_boot, .combine = "bind_rows") %dorng% { 
  df_boot <- boot_strat_newID(df_station, df_tmp, strat_n, "state", seed = i)
  df_out <- fit_break_discontinuous(station_year_df = df_boot, 
                                    fit_family_fn = use_fun, 
                                    # Unlike fit_break_trend, fitstat_dir is 
                                    # calculated within the function based on 
                                    # character string here
                                    break_fitstat = break_fitstat, 
                                    # Should be "mean" or "extremeN"
                                    PM_measure = pm_type, 
                                    break_alpha = break_alpha, 
                                    breaks_by = 1,  
                                    offset = NULL) %>%
    mutate(sample_name = sample_names, 
           boot = i)
  return(df_out)
}

write_feather(df_full, comb_name)

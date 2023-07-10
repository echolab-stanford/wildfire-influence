# ------------------------------------------------------------------------------
# Written by: Brandon de la Cuesta
# Conducts estimation using different hyperparameters.
# ------------------------------------------------------------------------------
# Set conus = T for US-wide breakpoint and F for regional
conus = T

# Set pm_type = "mean" for mean PM and "extremeN" (or anything not "mean") for extreme PM
pm_type = "extremeN"

# Set parallel = T if you want to use parallel processing (highly recommended), F otherwise
parallel = T

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

print(rep(paste0("CURRENTLY DOING: ", ifelse(conus, "CONUS", "REGIONAL"), "_", pm_type), 5))

# Setting overwrite = T will regenerate all subsamples
# Otherwise it will only do samples not present in directory
# Setting overwrite = F is useful when generating results and comparing where the 
# loop left off due to issues such as memory error
overwrite <- T

# Use cores_prop% of all cores
# Note: on Windows availableCores() will give you a number equal to the number 
# of threads for some CPUs
cores_prop = 0.35; n_cores <- round(cores_prop*availableCores())

# Condition check to set other globals/global fixest options
# If parallelizing, prevent fixest from using more than 1 thread at a time
if (parallel) {setFixest_nthreads(1)}
# If not parallelizing, use cores_prop% of cores
if (!parallel) {setFixest_nthreads(n_cores)}
if (pm_type == "extremeN") {offset <- as.formula("~log(n)")}
if (pm_type != "extremeN") {offset <- NULL}

# Overwrite globals as necessary
tol <- c(1e-5, 1e-5, 1e-8, 1e-8, 1e-10, 1e-10)
interval_width <- c(2, 4, 2, 4, 2, 4)
sample_obs <- rep(50, length(tol))
sample_years <- rep(15, length(tol))
sample_drop <- rep(0, length(tol))
sample_names <- c(paste0("tol_", tol, "_window", interval_width))
n_boot <- 500
file_names <- file.path(
  path_dropbox, "output", 
  paste0(ifelse(conus, "other/param_compare_conus_", 
                "other/param_compare_regional_"), 
         pm_type, "_", n_boot, "b", 
         "_alpha", break_alpha*100, "_corner", corner_push*100,  
         "_fitstat", break_fitstat, ifelse(smart_init, "_smartinit_", "_randominit_"), 
         "_maxit", max_it, "_maxrestart", max_restart, "_", sample_names, ".feather"))
comb_name <- file.path(
  path_dropbox, "output", 
  paste0(ifelse(conus, "other/param_compare_conus_", 
                "other/param_compare_regional_"), 
         pm_type, "_", n_boot, "b", 
         "_alpha", break_alpha*100,  "_corner", corner_push*100,  
         "_fitstat", break_fitstat, ifelse(smart_init, "_smartinit_", "_randominit_"), 
         "_maxit", max_it, "_maxrestart", max_restart, "_COMBINED.feather"))

# ------------------------------------------------------------------------------
#### Estimation ####
for (j in 1:length(sample_names)) {
  # Declare globals and backend
  if (parallel) {
    # Reset connections for each subsample to prevent dead workers
    closeAllConnections()
    registerDoFuture()
    plan(multisession, workers = n_cores)
  } 
  if (!parallel) {plan(sequential)}; tic()
  
  # Only run if the feather doesn't already exist
  if (!file.exists(file_names[j]) | overwrite) {
    
    # Subset to the sample you want
    # filter order is to restrict the sample based on desired station criteria 
    # before generating drop samples
    # This is done to simulate a case in which we use our primary station inclusion 
    # criteria but do not observe one of the later-period extreme years
    df_tmp <-  df_year %>% 
      filter(totalPM_n > sample_obs[j]) %>%
      group_by(id) %>% 
      mutate(n_years = n_distinct(year)) %>% 
      ungroup() %>%
      filter(n_years >= sample_years[j])  %>% 
      filter(year != sample_drop[j])
    
    # If doing regional, filter by the epa region as well
    if (!conus) {
      df_tmp <- df_tmp %>% filter(climate_regions == sample_names[j])
    }
    
    # Create station-level dataset, required for bootstrap function
    df_station <- df_tmp %>% 
      group_by(id, state) %>% 
      summarise() 
    
    # State-level station counts
    strat_n <- tapply(df_station$id, df_station$state, 
                      function(x){length(unique(x))})
    
    # Get smart initialization bounds
    break_init <- smart_initialize(
      station_year_df = df_tmp, 
      # Which fixest function you want to use
      # Current code set up to allow only feglm since it handles both ols and 
      # poisson in an MLE framework that tracks iterations and convergence
      fit_family_fn = use_fun, 
      # Current code registers "dev" as a fitstat
      break_fitstat = break_fitstat, 
      # Should be "mean" or "extremeN"
      PM_measure = pm_type, 
      # Choose to estimate piecewise or discontinuous
      # Always continuous for main results but retained for compatability in 
      # discontinuous script
      continuous_at_break = T, 
      # If you want to maximize set dir = -1, minimize dir = 1
      # Should always be 1 in this script since we're always using deviance as fitstat
      fitstat_dir = fitstat_dir, 
      # How much of the edges of the time variable you want to trim
      break_alpha = break_alpha, 
      # Default to half-year breaks
      # More breaks will give tighter smart initialization bounds 
      breaks_by = breaks_by, 
      # Time interval (width defined by breaks_by) you want to add to either 
      # side of the candidate breakpoint to build the uniform distribution
      interval_width = interval_width[j], 
      # Offset for poisson regression
      offset = offset, 
      # Link function
      family = ifelse(pm_type == "mean", "gaussian", "poisson"), 
      # Maximum iterations in fixest algorithm
      glm_it = glm_it
    )
    
    # Build list of (non-vectorized) arguments used in fit_break_trend and 
    # bootstrap function
    trend_args <- list(fit_family_fn = use_fun, 
                       PM_measure = pm_type, 
                       break_alpha = break_alpha, 
                       corner_push = corner_push, 
                       break_init = break_init, 
                       tol = tol[j], 
                       fe_fml = "id", 
                       max_it = max_it, 
                       max_restart = max_restart, 
                       offset = offset, 
                       glm_it = glm_it)
    boot_args <- list(df_ids = df_station, 
                      df_full = df_tmp, 
                      strat_n = strat_n, 
                      strat_var = "state")
    
    # Generate bootstrap samples and estimate
    list_boot <- future_Map(
      boot_strat_newID, 
      seed = 1:n_boot, 
      MoreArgs = boot_args, 
      # Turn off seed because bootstrap function calls this directly
      future.seed = NULL
    )
    list_sample <- future_Map(
      fit_break_trend, 
      station_year_df = list_boot, 
      # We need future.seed here because we're doing bootstrap 
      # and estimation separately, so need this for reproducibility 
      # of estimation results
      MoreArgs = trend_args, 
      future.seed = 1234, 
      # future.scheduling and future.chunk.size arguments can be set
      # to force single future for each element in X
      future.chunk.size = n_boot/(n_cores)
    )
    
    # Reshape to df and save
    df_sample <- do.call("bind_rows", list_sample)
    df_sample$sample_name <- sample_names[j]
    write_feather(df_sample, file_names[j])
    print(message("Finished with sample ", sample_names[j])); toc()
    
  } # Close if file.exists loop
}  # Close j loop

# Grab the files, bind, and save
df_full <- foreach(i = 1:length(file_names), .combine = "bind_rows") %do% {
  if (file.exists(file_names[i])) {
    return(read_feather(file_names[i]))
  }
}
write_feather(df_full, comb_name)

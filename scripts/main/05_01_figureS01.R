# ------------------------------------------------------------------------------
# Written by: Marissa Childs
# Plots Figure S1.
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

# Simulation results ----
param_comp <- feather::read_feather(file.path(path_dropbox, "output", "other", "param_compare_conus_mean_500b_alpha20_corner95_fitstatdev_smartinit__maxit5_maxrestart2_COMBINED.feather")) %>% 
  mutate(sample_name = ifelse(sample_name == "tol_1e-08_window2", 
                              paste0("Main (", sample_name, ")"), 
                              sample_name)) %>% 
  mutate(boot = 1:n(), .by = c(sample_name, state)) %>%
  # Filter out if either the total or nonsmoke regression didn't converge
  filter(if_all(starts_with("conv"), ~ .x) & if_all(matches("^b1|^b2"), ~!is.na(.x))) %>% 
  # Calculate coefficient differences
  mutate(b2nonsmokeMb2total = b2nonsmoke - b2total, 
         b2totalMb1total = b2total - b1total)

psi_test <- rbind(readRDS(file.path(path_dropbox, "output", "other", "algorithm_compare_gaussian.rds")) %>% 
                    mutate(family = "gaussian"), 
                  readRDS(file.path(path_dropbox, "output", "other", "algorithm_compare_poisson.rds")) %>% 
                    mutate(family = "poisson"))

# ------------------------------------------------------------------------------
# Figure S1: algorithm ------
plot_grid(
  plot_grid(
    # Panel a: results in simulation
    psi_test %>% 
      ggplot(aes(x = psi_true, y = psi, color = test, shape = test)) + 
      geom_point() + 
      geom_abline(slope = 1, intercept = 0) + 
      facet_wrap(~family, 
                 nrow = 2) + 
      scale_color_manual(name = "", 
                         values = c("red", "blue"), 
                         labels = c("ours", "segmented")) + 
      scale_shape_manual(name = "", 
                         values = c(16, 5), 
                         labels = c("ours", "segmented")) + 
      xlab("true breakpoint") + ylab("estimated breakpoint") + 
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = c(0.18, 0.9), 
            legend.title=element_blank(),
            legend.background = element_blank(), 
            legend.box.margin=margin(-10,0,-4,-5),
            legend.box.background = element_rect(color = "grey10")),
    
    # Panel b: distribution of breaks by hyperparams 
    param_comp %>% 
      dplyr::select(break_year, boot, sample_name) %>% 
      separate(sample_name, into = c(NA, "tol", "window"), sep = "tol_|_window", 
               remove = FALSE) %>% 
      mutate(window = paste0("initialization window: ", gsub(")", "", window)), 
             tol = paste0("convergence\ntol: ", tol)) %>% 
      unique %>% 
      ggplot(aes(x = break_year)) + 
      geom_histogram(aes(y = after_stat(density))) + 
      theme_classic() + 
      facet_grid(window ~ tol) + 
      theme(strip.background = element_blank(), 
            text = element_text(size = 8), 
            strip.text = element_text(size = 10)) + 
      xlab("breakpoint year"), 
    nrow = 2, 
    align = "v",
    axis = "l",
    labels = c("a", "b"), 
    rel_heights = c(1, 0.9)),
  
  # Panel c: classifications by hyperparameters
  param_comp %>% 
    dplyr::select(state, sample_name, break_year, starts_with("b"), -boot) %>% 
    group_by(state, sample_name) %>%
    classify_boot(full = T) %>% 
    mutate(sample_name = gsub("tol_", "tol ", sample_name) %>% 
             gsub("_window", ", window ", .)) %>% 
    {ggplot(data = ., 
            aes(y = state_name,
                x = sample_name, 
                color = smoke_group_factor, 
                fill = smoke_group_factor)) + 
        geom_tile() + 
        scale_y_discrete(limits = rev) + 
        scale_color_manual("Smoke-influence classification",
                           values = smokeinfluence_colors,
                           aesthetics = c("color", "fill")) +
        ylab("") + xlab("") + 
        theme_classic() + 
        theme(axis.text = element_text(size = 8),
              axis.text.x = element_text(angle = 30, hjust = 1),
              legend.position = "bottom", 
              legend.direction = "vertical", 
              legend.key.size = unit(15, "points"),
              legend.text = element_text(size = 8),
              legend.box.margin=margin(-37,-10,0,-10),
              plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points"))}, 
  nrow = 1, 
  rel_widths = c(1, 0.8),
  labels = c("", "c")) %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureS01_algorithm_performance.pdf"), ., 
         width = 8, height = 9) 

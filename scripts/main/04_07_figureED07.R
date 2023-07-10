# ------------------------------------------------------------------------------
# Written by: Marissa Childs
# Plots Figure ED7.
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

# ------------------------------------------------------------------------------
# Figure ED7: extreme days sensitivity ----
plot_grid(
  # Panel a
  extreme_classified %>% 
    filter(!(state_abbr == "US" & grepl("region", sample_name))) %>% 
    {ggplot(data = ., 
            aes(y = state_name,
                x = sample_factor, 
                color = smoke_influenced_text, 
                fill = smoke_influenced_text)) + 
        geom_tile() + 
        scale_y_discrete(limits = rev) + 
        scale_color_manual("Smoke-influence classification",
                           values = c("grey90", "#5a5a83"),
                           aesthetics = c("color", "fill")) +
        ylab("") + xlab("") + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1), 
              legend.position = "none")}, 
  # Panel b
  extreme_classified %>% 
    filter(state_abbr != "US") %>% 
    {ggplot(data = ., 
            aes(y = sample_factor, 
                color = smoke_influenced_text, 
                fill = smoke_influenced_text)) + 
        geom_bar() + 
        scale_color_manual("Smoke-influence classification",
                           values = extreme_colors,
                           aesthetics = c("color", "fill")) + 
        ggnewscale::new_scale_color() + 
        geom_text(data = extreme_classified %>% 
                    filter(state_abbr != "US") %>% 
                    summarise(n = n(), .by = c(sample_factor, smoke_influenced_text)) %>% 
                    arrange(sample_factor, desc(smoke_influenced_text)) %>% 
                    mutate(lower = cumsum(replace_na(lag(n), 0)),
                           upper = cumsum(n),
                           mid = (lower + upper)/2,
                           .by = sample_factor), 
                  aes(x = mid, y = sample_factor, label = n, 
                      color = I(ifelse(smoke_influenced_text == "smoke-influenced", "grey90", "grey10"))), 
                  inherit.aes = F) + 
        scale_y_discrete(limits = rev) +
        xlab("Number of states") + ylab("") + 
        scale_x_continuous(breaks = seq(0, 48, by = 8)) + 
        theme_classic() + 
        theme(legend.position = "bottom", 
              legend.direction = "vertical")}, 
  ncol = 2,
  rel_widths = c(0.8, 1),
  labels = "auto") %>%
  ggsave(file.path(path_github, "figures", "raw","figureED07_extreme_smokeinfluenced_sensitivity.pdf"), ., 
         width = 8, height = 7)

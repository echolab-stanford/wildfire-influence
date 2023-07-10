# ------------------------------------------------------------------------------
# Written by: Marissa Childs
# Plots Figure ED8.
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
                       # add row to get down to zero 
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

# ------------------------------------------------------------------------------
# Figure ED8: extreme days from wildfire smoke by year/period ----
plot_grid(
  # Panel a
  extremePct_cdf %>% 
    filter(type == "annual") %>%
    ggplot(aes(x = pct_exceed,
               y = period, 
               group = period, 
               color = period,
               fill = period)) +
    geom_density_ridges(alpha = 0.6, color = alpha("grey20", 0.7),
                        rel_min_height = 0.01,
                        bandwidth = 0.04,
                        jittered_points = TRUE,
                        position = position_points_jitter(width = 0.005, 
                                                          height = 0),
                        point_color = "black", point_alpha = 0.5,
                        point_shape = "|", point_size = 1.5) +
    scale_color_manual(name = "year", 
                       values = rev(MetBrewer::met.brewer("Hiroshige", 17)), 
                       aesthetics = c("fill", "color")) +
    xlab(expression(paste("% of days > 35 ", mu, "g/", m^3, " due to smoke"))) + 
    scale_y_discrete(expand = expansion(mult = c(0.03, 0.075))) +
    scale_x_continuous(lim = c(-0.001, 1), 
                       expand = expansion(mult = c(0.02, 0.05))) +
    ylab("") + 
    theme_classic() + 
    theme(legend.position = "none", 
          plot.margin = unit(c(5.5, 5.5, 5.5, 2.5), "pt")), 
  # Panel b 
  extremePct_cdf %>%   
    filter(type != "annual" | period == "2021") %>% 
    {ggplot(mapping = aes(x = pct_exceed, y = n_gte, 
                          group = period)) + 
        geom_step(data = filter(., group == "annual"),
                  mapping = aes(color = period),
                  lwd = 0.6, direction = "vh") +
        geom_step(data = filter(., group != "annual"),
                  mapping = aes(color = period),
                  lwd = 1, direction = "vh") +
        scale_color_manual(name = "",
                           values = c(rev(MetBrewer::met.brewer("Hiroshige", 7)[c(1, 3, 7)]), "grey30")) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
        scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)) ,
                           labels = scales::percent) +
        xlab(expression(paste("% of days > 35 ", mu, "g/", m^3, " due to smoke"))) + 
        ylab("Number of states") + 
        theme_classic() + 
        theme(legend.key.height = unit(13, "pt"), 
              legend.box = "horizontal", 
              legend.position = c(0.65, 0.85), 
              legend.justification = c(0, 0.5), 
              plot.margin = unit(c(5.5, 8.5, 5.5, 5.5), "pt"))}, 
  labels = "auto", rel_widths = c(0.8, 1)) %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureED08_extremePct_from_smoke.pdf"), ., 
         width = 7, height = 4.5)

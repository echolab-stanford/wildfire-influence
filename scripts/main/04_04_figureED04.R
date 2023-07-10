# ------------------------------------------------------------------------------
# Written by: Marissa Childs
# Plots Figure ED4.
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
# Figure ED4: coefficient differences sensitivity ----
bind_rows(mean_classified %>% 
            dplyr::select(state_name, sample_name, sample_factor, contains("Mb")) %>% 
            mutate(type = "mean"),
          extreme_classified %>% 
            dplyr::select(state_name, sample_name, sample_factor, contains("Mb")) %>% 
            mutate(type = "extreme")) %>% 
  filter(!(state_name == "Contiguous US" & grepl("region", sample_name))) %>% 
  pivot_longer(starts_with("b"), 
               names_sep = "_", 
               names_to = c("par", ".value")) %>% 
  filter(par == "b2nonsmokeMb2total" | type == "mean") %>% 
  mutate(spec_type = case_when(sample_name == main_spec ~ "main",
                               grepl('drop', sample_name) ~ "dropped years",
                               grepl("^n|all", sample_name) ~ "sample restriction",
                               grepl("regional", sample_name) ~ "region breaks",
                               sample_name == "discontinuous" ~ "discontinuous", 
                               sample_name == "pos_anom" ~ "positive smoke\nPM2.5 anomalies",
                               T ~ "other"),
         spec_type = factor(spec_type,
                            levels = c("other",
                                       "positive smoke\nPM2.5 anomalies",
                                       "discontinuous", 
                                       "region breaks",
                                       "dropped years",
                                       "sample restriction",
                                       "main"),
                            ordered = T)) %>% 
  left_join(data.frame(par = c("b2totalMb1total", "b2nonsmokeMb2total")) %>% 
              mutate(par_exp = gsub("M", " - ", par) %>% 
                       gsub("total", "", .) %>% 
                       gsub("nonsmoke", "\'", .))) %>% 
  mutate(panel = paste0(type, ", ", par_exp) %>% 
           factor(levels = sort(unique(.))[c(2,3,1)], ordered = T)) %>% 
  {ggplot(data = ., 
          aes(y = as.factor(state_name), 
              x = mid, 
              color = spec_type,
              group = spec_type)) +
      geom_vline(xintercept = 0) +
      geom_errorbar(aes(xmin = lower,
                        xmax = upper),
                    width = 0,
                    position = position_dodge(width = 0.6)) +
      geom_point(aes(size = I(ifelse(sample_name == main_spec, 2, 1.5)),
                     shape = I(ifelse(sample_name == main_spec, 1, 16))),
                 position = position_dodge(width = 0.6)) + 
      facet_wrap(~panel, scales = "free_x", 
                 labeller = as_labeller(
                   c(`mean, b2 - b1` = "\u0394 \u03bc g/m\u00B3/year", 
                     `mean, b2' - b2` = "\u0394 \u03bc g/m\u00B3/year", 
                     `extreme, b2' - b2` = "% \u0394 in pct \nextreme days/year")), 
                 strip.position = "bottom") + 
      scale_color_manual(values = c("#a82203", "#208cc9","purple3", "#f1af3a",
                                    "#637b31", "black"),
                         name = "specification type",
                         guide = guide_legend(reverse = TRUE,
                                              override.aes = list(shape = c(1, 16, 16, 16, 16, 16)))) +
      scale_y_discrete(limits = rev) +
      coord_cartesian(clip = "off") +
      theme_classic() + 
      theme(strip.background = element_blank(), 
            strip.placement = "outside", 
            strip.text = element_text(size = 12), 
            panel.spacing.x = unit(16.5, "pt"), 
            plot.margin = unit(c(32.5, 5.5, 5.5, 5.5), "points")) + 
      ylab("") + 
      xlab("")}  %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureED04_coefdiff_sensitivity.pdf"), ., 
         width = 7.5, height = 7, 
         device = cairo_pdf)

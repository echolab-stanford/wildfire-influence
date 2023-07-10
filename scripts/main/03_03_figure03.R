# ------------------------------------------------------------------------------
# Written by: Marissa Childs
# Plots Figure 3.
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
# Figure 3: mean small multiples ---- 
purrr::pmap(
  data.frame(plot_spec = main_spec,
             class_df = c("mean_classified", "extreme_classified"),
             pm_measure = c("mean", "extremePct"),
             grouping = c("smoke_group_factor", "smoke_influenced_text"),
             class_colors = I(list(alpha(smokeinfluence_colors, 
                                         c(0.6, 0.7, 0.8, 0.8, 0.8, 0.4)), 
                                   alpha(extreme_colors, c(0.4, 0.6)))), 
             ylab_fn = I(c(NA, scales::percent)), 
             fig_name = c("figure03_smoke_influence.pdf", 
                          "figureED03_extremePct_smoke_influence.pdf"))[1,],
  function(plot_spec, class_df, pm_measure, grouping, class_colors, 
           ylab_fn, fig_name){
    if(is.na(ylab_fn)){
      y_lab_fn <- waiver()
    } else{
      y_lab_fn <- ylab_fn
    }
    fig <- state_year %>% 
      rename(value = all_of(pm_measure)) %>%
      ggplot() + 
      geom_rect(data = get(class_df) %>%
                  filter(sample_name == plot_spec) %>%
                  rename(color_by = all_of(grouping)),
                aes(fill = color_by),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf) +
      geom_line(aes(x = year, y = value,
                    group = interaction(state_abbr, pm_type), 
                    color = pm_type)) + 
      geom_vline(data = get(class_df) %>%
                   filter(sample_name == plot_spec),
                 aes(xintercept = break_year_mid), 
                 linetype = "dotted") +
      facet_geo(~state_abbr, scales = "free_y",
                grid = my_grid, label = "code") +
      scale_color_manual(values = c("blue",
                                    "grey10"),
                         guide = "none") +
      scale_fill_manual(name = "",
                        values = class_colors) +
      theme_classic(base_size = 9.5) +
      theme(strip.background = element_rect(linetype = "blank", fill = NA),
            strip.text = element_text(face = "bold",
                                      size = 10.5,
                                      vjust = -0.5),
            panel.spacing.x = unit(3, "points"),
            panel.spacing.y = unit(-2, "points"),
            legend.position = c(0.87, 0.2),
            legend.text = element_text(size = 11),
            axis.title = element_text(size = 12.5),
            plot.margin = unit(c(5.5, 8, 4, 5.5), "pt")) +
      scale_x_continuous(breaks = c(2000, 2022),
                         labels = c("   2000", "2022   ")) +
      scale_y_continuous(breaks = scales::breaks_extended(n = 4, Q = c(1, 2, 3, 4, 5, 6)), 
                         labels = y_lab_fn) +
      xlab("") + 
      {if(pm_measure == "mean"){ylab(c(expression(paste(PM[2.5]," (",mu, "g/", m^3, ")"))))
      }else if(pm_measure == "extremePct"){ylab(expression(paste("% of days > 35 ",mu, "g/", m^3)))}}
    
    # Make it a gtable
    fig %<>% ggplot_build %>% ggplot_gtable
    
    fig %<>% clean_geo_facet()
    
    ggsave(file.path(path_github, "figures", "raw", fig_name), 
           fig, width = 11.5, height = 8.5)
    
  })

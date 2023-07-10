# ------------------------------------------------------------------------------
# Written by: Marissa Childs
# Plots Figure 2.
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

# Calculate percent undone ----
pct_undone <- mean_fits %>%
  filter(state != "US" & sample_name == main_spec) %>% 
  mutate(boot = 1:n(), 
         .by = c(sample_name, state)) %>% 
  left_join(mean_classified %>% 
              filter(sample_name == main_spec) %>% 
              dplyr::select(state = state_name, 
                            smoke_influenced, smoke_group)) %>% 
  filter(smoke_influenced) %>% dplyr::select(-smoke_influenced) %>% 
  transmute(boot, smoke_group, state,
            changeEarlyTotal = b1total*(break_year - 2000),
            changeLateSmoke = (-b2nonsmokeMb2total)*(2022 - break_year),
            pctUndone = pmax(pmin(-changeLateSmoke/changeEarlyTotal, 1), 0)*100, 
            yearsUndone = changeLateSmoke/b1total)
pct_undone_groups <- pct_undone %>% 
  dplyr::select(-state) %>% 
  {rbind(., 
         mutate(., smoke_group = "all_smoke-influenced"), 
         filter(., grepl("reversal", smoke_group)) %>% 
           mutate(smoke_group = "smoke_reversals"))} %>% 
  # For each bootstrap, calculate avg (median) pct undone, 
  summarise(across(everything(), 
                   list(mid =~ median(.x), 
                        mean =~ mean(.x))), 
            .by = c(boot, smoke_group)) %>% 
  summarise(across(!boot, 
                   list(lower =~ quantile(.x, q_alpha/2), 
                        mid =~ quantile(.x, 0.5),
                        # mean =~ quantile(.x, 0.5),
                        upper =~ quantile(.x, 1 - q_alpha/2))), 
            .by = smoke_group) %>% 
  pivot_longer(!smoke_group) %>% 
  separate_wider_delim(cols = name, 
                       names = c("name", "sample_stat", "boot_stat"), 
                       delim = "_") %>% 
  mutate(panel = ifelse(grepl("^change", name), 
                        "change in PM2.5",
                        "pct of progress undone by smoke"))

# ------------------------------------------------------------------------------
# Figure 2a: sankey of classifications ----
mean_classified %>% 
  filter(sample_name == main_spec & state_abbr != "US") %>% 
  group_by(smoke_group, change_group, smoke_group_factor, change_group_factor) %>% 
  summarise(count = n(), .groups = "drop") %>%  
  gather_set_data(1:2) %>% 
  {mutate(., across(c(y, smoke_group, change_group),
                    function(x) factor(x,
                                       levels = c(levels(.$change_group_factor), levels(.$smoke_group_factor)),
                                       ordered = TRUE)))} %>% 
  dplyr::select(-ends_with("factor")) %>% 
  group_by(y) %>%
  mutate(group_count = sum(count)) %>%
  ungroup %>%
  {ggplot(data = ., 
          aes(y = x, id = id, split = y, value = count)) +
      geom_parallel_sets(aes(fill = smoke_group), 
                         alpha = 0.3, axis.width = 0.1) +
      geom_parallel_sets_axes(aes(fill = y), 
                              axis.width = 0.1) +
      geom_text(data = mutate(., lab_position_y = ifelse(y == change_group, 2, 1)) %>% 
                  dplyr::select(y, group_count, lab_position_y) %>% 
                  unique %>%
                  group_by(lab_position_y) %>% 
                  arrange(desc(y)) %>% 
                  mutate(lab_position_x = 2.5*((1:n()) - 1) + cumsum(group_count) - group_count/2),
                aes(x = lab_position_x, 
                    y = lab_position_y, 
                    label = group_count), 
                size = 3,
                inherit.aes = FALSE) +
      geom_text(data = mutate(., lab_position_y = ifelse(y == change_group, 2, 1)) %>% 
                  dplyr::select(y, group_count, lab_position_y) %>% 
                  unique %>%
                  group_by(lab_position_y) %>% 
                  arrange(desc(y)) %>% 
                  mutate(lab_position_x = 2.5*((1:n()) - 1) + cumsum(group_count) - group_count/2) %>% 
                  mutate(axis_lab = gsub("^smoke-", "", y) %>% 
                           gsub(" stagnation", "\nstagnation", .) %>% 
                           gsub(" reversal", "\nreversal", .) %>% 
                           gsub("acceleration", "accel.", .) %>%  
                           gsub("influence ", "influence\n",. )  %>% 
                           gsub("early decline", "early\ndecline",. ) %>% 
                           gsub("^non-sig.", "non-sig.\n", .) %>% 
                           gsub("trend change", "trend\nchange", .)) %>% 
                  filter(lab_position_y == 1),
                aes(x = lab_position_x, 
                    y = lab_position_y, 
                    label = axis_lab), 
                vjust = 1, 
                size = 2.6,
                fontface = "italic",
                lineheight = 0.8,
                nudge_y = -0.08,
                inherit.aes = FALSE) +
      geom_text(data = mutate(., lab_position_y = ifelse(y == change_group, 2, 1)) %>% 
                  dplyr::select(y, group_count, lab_position_y) %>% 
                  unique %>%
                  group_by(lab_position_y) %>% 
                  arrange(desc(y)) %>% 
                  mutate(lab_position_x = 2.5*((1:n()) - 1) + cumsum(group_count) - group_count/2) %>% 
                  mutate(axis_lab = gsub("^smoke-", "", y) %>% 
                           gsub(" stagnation", "\nstagnation", .) %>% 
                           gsub(" reversal", "\nreversal", .) %>% 
                           gsub("acceleration", "accel.", .) %>%  
                           gsub("influence ", "influence\n",. )  %>% 
                           gsub("early decline", "early\ndecline",. ) %>% 
                           gsub("^non-sig.", "non-sig.\n", .) %>% 
                           gsub("trend change", "trend\nchange", .)) %>% 
                  filter(lab_position_y == 2),
                aes(x = lab_position_x, 
                    y = lab_position_y, 
                    label = axis_lab), 
                vjust = 0, 
                fontface = "italic",
                nudge_y = 0.08,
                size = 2.6,
                inherit.aes = FALSE) +
      scale_fill_manual(values = c("#c38f16",
                                   "#5d6174", # non-sig change
                                   "grey60", # not smoke influenced
                                   "#dcdcdb", # no sig early decline 
                                   alpha("#175449", 0.85), # reversal
                                   alpha("#62205f", 0.8), # smoke-caused reversal
                                   "#bb292c", # smoke-influenced reversal
                                   "#ef8737", # smoke-influenced stagnation
                                   "#ffd353", # smoke-influenced, no early decline
                                   "#dea868",
                                   "#6e948c")) +  # stagnation 
      scale_y_continuous(expand = expansion(mult = 0.2)) +
      xlab("") + ylab("") + theme_void() + 
      theme(legend.position = "none")} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figure02a_classification_sankey.pdf"), ., 
         width = 6.5, height = 2.5)

# Figure 2b: pct undone ---- 
pct_undone %>% 
  dplyr::select(-smoke_group, -boot, -yearsUndone) %>% 
  summarise(across(everything(), median), 
            .by = state) %>% 
  pivot_longer(!state) %>% 
  mutate(panel = ifelse(grepl("^change", name), 
                        "change in PM2.5",
                        "pct of progress undone by smoke")) %>% 
  {ggplot(data = ., 
          aes(x = value, 
              group = name,
              fill = name)) + 
      geom_histogram(alpha = 0.7, 
                     color = "black",
                     position = "identity") + 
      geom_vline(data = pct_undone_groups %>% 
                   filter(smoke_group == "all_smoke-influenced") %>% 
                   filter(boot_stat %in% c("mid", "mean") & 
                            name != "yearsUndone"),
                 aes(xintercept = value,
                     linetype = sample_stat,
                     color = name)) + 
      facet_wrap(~panel, 
                 scales = "free", nrow = 2, 
                 strip.position = "bottom") + 
      scale_fill_manual(values = c("lightblue", "orange", "grey")) +
      scale_color_manual(values = c("lightblue", "black", "red")) +
      scale_linetype_manual(values = c("solid", "dashed")) + 
      ylab("# of states") + 
      xlab("") + 
      theme_classic() + 
      theme(strip.background = element_blank(), 
            strip.placement = "outside", 
            strip.text = element_text(size = 11),
            legend.position = "none",
            axis.title = element_text(size = 11))} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figure02b-c_pct_undone.pdf"), ., 
         width = 4, height = 3.5)

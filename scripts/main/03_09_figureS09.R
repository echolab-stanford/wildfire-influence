source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

# set up state grid
my_grid <- us_state_grid3
# move ME up, VT and NH one right, FL one left
my_grid %<>% 
  mutate(row = ifelse(code == "ME", row - 1, row), 
         col = ifelse(code %in% c("VT", "NH"), col + 1, col), 
         col = ifelse(code == "FL", col - 1, col)) %>% 
  # remove AK and HI
  filter(code %in% c("AK", "HI") == FALSE) %>% 
  # add US to top left
  rbind(data.frame(row = 1, col = 1, code = "US", name = "Contiguous US"))

# state dictionary ----
epa_regions <- read.csv(file.path(path_dropbox, "data", "us_climate_regions.csv"))
nonContig_stateFIPS <- c("02","60","66","15","72","78","69")

state_dict <- tigris::states(cb =  TRUE) %>% 
  filter(STATEFP %in% nonContig_stateFIPS == FALSE) %>% 
  left_join(epa_regions, by = c("STUSPS" = "state_code")) %>% 
  mutate(climate_regions = gsub(" \\(.*", "", climate_regions),
         climate_regions = stringr::str_to_title(climate_regions)) %>% 
  sf::st_drop_geometry() %>% 
  select(state_code = STATEFP, 
         state_abbr = STUSPS, 
         state_name = NAME,
         climate_regions = climate_regions)

# station-year data ----
station_year <- readRDS(file.path(path_dropbox, "data", 
                                  "epa_station_year_filled2022.rds")) %>% 
  rename(state_name = state) %>% 
  left_join(state_dict)

# state-year averages ----
state_year <- station_year %>% 
  filter(totalPM_n > 50) %>% 
  group_by(id) %>% 
  mutate(n_year = n()) %>% 
  ungroup %>% 
  filter(n_year >= 15) %>% 
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

# read in fits ----
all_fits <- readRDS(file.path(path_dropbox, "output", 
                              "state_allFits_coefs_classifications.rds")) %>% 
  mutate(smoke_group = gsub("no smoke-influence detected", 
                            "no smoke influence detected", 
                            smoke_group)) %>% 
  mutate(smoke_group = gsub(", ", ",\n", smoke_group),
         smoke_group_factor = factor(smoke_group, 
                                     levels = c("smoke-caused reversal", 
                                                "smoke-influenced reversal", 
                                                "smoke-influenced stagnation", 
                                                "smoke-influenced,\nno early decline", 
                                                "smoke-influenced,\nno signficant trend change", 
                                                "no smoke influence detected"), 
                                     ordered = TRUE), 
         change_group = gsub("sig", "sig.", change_group),
         change_group_factor = factor(change_group, 
                                      levels = c("reversal", 
                                                 "stagnation",
                                                 "acceleration", 
                                                 "non-sig. change", 
                                                 "no sig. early decline"), 
                                      ordered = TRUE),
         smoke_influenced_text = ifelse(smoke_influenced, 
                                        "smoke-influenced", 
                                        "no smoke influence detected"), 
         spec_factor = ifelse(spec %in% c("break_1early", "break_1late") & 
                                PM_measure == "extremePct", 
                              paste0(spec, "_", PM_measure), 
                              spec), 
         spec_factor = recode_factor(spec_factor, 
                                     "main" = "Main", 
                                     "allobs" = "All obs", 
                                     "50obs5yr" = "5yrs, 50obs",
                                     "50obs10yr" = "10yrs, 50obs",
                                     "50obs15yr" = "15yrs, 50obs", 
                                     "100obs5yr" = "5yrs, 100obs", 
                                     "100obs10yr" = "10yrs, 100obs", 
                                     "break_1early" = "break in 2015", 
                                     "break_1late" = "break in 2017", 
                                     "break_1early_extremePct" = "break in 2009", 
                                     "break_1late_extremePct" = "break in 2011", 
                                     "drop2020" = "drop 2020", 
                                     "drop2021" = "drop 2021", 
                                     "drop2022" = "drop 2022", 
                                     "no_duplicated_break" = "no dup. break yr", 
                                     "piecewise" = "piecewise", 
                                     "regional_breaks" = "regional breaks",
                                     .ordered = TRUE))

# fig S9: western states sensitivity to break years -----
breaks <- rbind(read.csv(file.path(path_dropbox, "output", "year_subsample_estimates.csv")) %>%  
                  mutate(PM_measure = "mean"), 
                read.csv(file.path(path_dropbox, "output", "epm_year_subsample_estimates.csv")) %>% 
                  mutate(PM_measure = "extremePct")) %>%
  select(-X, -climate_regions) %>% 
  left_join(state_dict %>% 
              mutate(state_code = as.numeric(state_code))) %>% 
  select(-starts_with("state")) %>% 
  unique %>%
  rename_with(function(x) paste0("breakpoint-", x), 
              !c(climate_regions, PM_measure, starts_with("nrow"))) %>% 
  rename_with(function(x) gsub("nrow_", "nrow-", x), 
              starts_with("nrow")) %>% 
  pivot_longer(!c(climate_regions, PM_measure), 
               names_sep = "-", 
               names_to = c(".value", "spec"))

west_fits <- readRDS(file.path(path_dropbox, "output", 
                               "west_breakpointFits_coefs_classifications.rds"))

west_fits %>% 
  # mutate(trend = ifelse(b2totalMb1total_est > 0 & b2totalMb1total_pval < 0.05, "yes", "no"), 
  #        smoke = ifelse(b2nonsmokeMb2total_est < 0 & b2totalMb1total_pval < 0.05, "yes", "no")) %>%  
  pivot_longer(starts_with("b"), 
               names_sep = "_", 
               names_to = c("par", ".value")) %>% 
  mutate(par_clean = gsub("total", "", par) %>% 
           gsub("nonsmoke", "'", .)) %>%
  # filter(state_name %in% c("Washington", "Idaho", "Oregon", "Nevada", "California")) %>% 
  mutate(state_name = factor(state_name, 
                             levels = c("Washington", "Idaho", "Montana", 
                                        "Oregon", "Nevada", "Wyoming", 
                                        "California", "Utah", "Colorado"),
                             ordered = TRUE)) %>% 
  mutate(est = ifelse(par == "b2nonsmokeMb2total", -est, est)) %>%
  {ggplot() + 
      geom_vline(xintercept = 0) + 
      geom_errorbar(data = select(., state_name, climate_regions) %>% 
                      unique %>% 
                      left_join(breaks %>% 
                                  filter(spec == "n50_y10" & PM_measure == "mean")), 
                    aes(y = breakpoint, 
                        xmin = -Inf, xmax = 0.35), 
                    color = "black",
                    linetype = "dashed") + 
      geom_errorbar(data = filter(., par %in% c("b1total", "b2total", "b2nonsmoke")),
                    mapping = aes(y = split,
                                  xmin = est + qt(0.025, t_df)*se,
                                  xmax = est + qt(0.975, t_df)*se,
                                  color = par_clean),
                    width = 0) +
      geom_point(data = filter(., par %in% c("b1total", "b2total", "b2nonsmoke")),
                 mapping = aes(y = split, x = est, color = par_clean)) + 
      geom_text(data = filter(., grepl("Mb", par)),
                aes(x = ifelse(par == "b2totalMb1total", 0.47, 0.7), 
                    y = split,
                    label = ifelse(est > 0 & pval < 0.05, "yes", "no"))) +
      geom_text(data = data.frame(x = c(0.47, 0.7), 
                                  lab = c("b1<b2","b2'<b2")), 
                mapping = aes(x = x, label = lab), 
                y = 2017, size = 3, fontface = "bold") +
      # geom_text(data = ., 
      #           aes(x = 0.7, y = split, 
      #               label = ifelse(b2nonsmokeMb2total_est < 0 & b2totalMb1total_pval < 0.05, "yes", "no"))) + 
      # annotate("text", x = 0.7, y = 2017, label = "b2'<b2", 
      #          size = 3, fontface = "bold") +
      facet_wrap(~state_name, ncol = 2) + 
      coord_cartesian(clip = "off") + 
      scale_y_continuous(breaks = seq(2008, 2016, by = 2), 
                         expand = expansion(mult = c(0.1, 0.05))) + 
      scale_color_manual("", values = c("blue", "red", "orange")) + 
      # scale_x_continuous(limits = c(-0.5, NA)) + 
      xlab(expression(paste(PM[2.5]," trend (",mu, "g/", m^3, "/year)"))) +
      ylab("") + 
      theme_classic() + 
      theme(strip.background = element_blank(), 
            panel.spacing.x = unit(15, "pt"),
            plot.margin = unit(c(5.5, 8.5, 5.5, 3.5), "pt"), 
            legend.position = c(0.8, 0.15), 
            legend.background = element_rect(fill = NA))} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureS09_west_coef_breaksensitivity.pdf"), 
         ., width = 6, height = 7)

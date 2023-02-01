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

# fig S6: coefficient differences sensitivity ----

par_exp_df <- data.frame(par = c("b2totalMb1total", "b2nonsmokeMb2total")) %>% 
  mutate(par_exp = gsub("M", " - ", par) %>% 
           gsub("total", "", .) %>% 
           gsub("nonsmoke", "\'", .))
# mutate(par_exp = factor(par, 
#                         levels = c("b2totalMb1total", "b2nonsmokeMb2total"), 
#                         ordered = TRUE, 
#                         labels = c(expression(paste("b2 - b1 (", Delta, " ",mu, "g/", m^3, "/year)")), 
#                                    expression(paste("b2' - b2 (", Delta, " ",mu, "g/", m^3, "/year)")))))

all_fits %>%
  filter(PM_measure == "mean") %>%
  select(state_name, t_df, spec, starts_with("b"), -breakpoint) %>% 
  pivot_longer(starts_with("b"), 
               names_sep = "_", 
               names_to = c("par", ".value")) %>% 
  filter(grepl("Mb", par)) %>%
  mutate(spec_type = case_when(spec == "main" ~ "main", 
                               spec %in% c("no_duplicated_break", "piecewise") ~ "fitting", 
                               grepl("break_1", spec) ~ "break year", 
                               grepl('drop', spec) ~ "dropped years", 
                               grepl("obs|yr", spec) ~ "sample restriction", 
                               spec == "regional_breaks" ~ "region breaks",
                               T ~ "other"), 
         spec_type = factor(spec_type, 
                            levels = c("fitting", 
                                       "region breaks", 
                                       "break year", 
                                       "dropped years", 
                                       "sample restriction", 
                                       "main"), 
                            ordered = T)) %>%  
  left_join(par_exp_df) %>%
  arrange(spec_type) %>% 
  {ggplot(data = ., 
          aes(y = as.factor(state_name), 
              x = est, 
              color = spec_type, 
              group = spec)) + 
      geom_vline(xintercept = 0) +
      geom_errorbar(aes(xmin = est + qt(0.025, t_df)*se,
                        xmax = est + qt(0.975, t_df)*se),
                    width = 0, 
                    position = position_dodge(width = 0.6)) +
      geom_point(aes(size = I(ifelse(spec == "main", 2, 1.5)),
                     shape = I(ifelse(spec == "main", 1, 16))),
                 position = position_dodge(width = 0.6)) + 
      facet_grid(~par_exp, scales = "free") + 
      # labeller = label_parsed, switch = "x") + 
      # labeller = label_bquote(cols = Delta^.(par))) + 
      # labeller = list("b2totalMb1total" = expression(paste("b2 - b1 (", delta, ")")),
      # "b2nonsmokeMb2total" = expression(paste("b2' - b2 (", delta, ")")))) + 
      scale_color_manual(values = c("#a82203", "#208cc0", "#f1af3a",
                                    "#cf5e4e", "#637b31", "black"),
                         name = "specification type",
                         guide = guide_legend(reverse = TRUE, 
                                              override.aes = list(shape = c(1, 16, 16, 16, 16, 16)))) +
      scale_y_discrete(limits = rev) +
      theme_classic() + 
      theme(strip.background = element_blank(), 
            strip.placement = "outside", 
            strip.text = element_text(size = 12, face = "bold"), 
            panel.spacing.x = unit(16.5, "pt")) + 
      ylab("") + 
      xlab(expression(paste(Delta, " ",mu, "g/", m^3, "/year")))}  %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureS06_coefdiff_sensitivity.pdf"), 
         ., width = 7, height = 7)

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

# fig S7: trend-change classes sensitivity ----
pmtrend_colors <- c("#175449", "#6e948c", "#c38f16", "#dcdcdb", 
                    # alpha("#b9b9b8", 0.6), 
                    "#5d6174")
plot_grid(
  # panel a
  all_fits %>% 
    filter(PM_measure == "mean") %>% 
    select(state_name, spec_factor, change_group_factor) %>% 
    # rbind(data.frame(state_name = "Contiguous US", 
    #                  spec_factor = "regional breaks", 
    #                  change_group_factor = NA)) %>%
    {ggplot(data = ., 
            aes(y = state_name,
                x = spec_factor, 
                color = change_group_factor, 
                fill = change_group_factor)) + 
        geom_tile() + 
        scale_y_discrete(limits = rev) + 
        scale_color_manual(values = pmtrend_colors,
                           na.value = "white",
                           aesthetics = c("color", "fill")) +
        ylab("") + xlab("") + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1), 
              legend.position = "none")}, 
  # panel b
  all_fits %>% 
    filter(PM_measure == "mean") %>% 
    filter(state_abbr != "US") %>% 
    {ggplot(data = ., 
            aes(y = spec_factor, 
                color = change_group_factor, 
                fill = change_group_factor)) + 
        geom_bar() + 
        scale_y_discrete(limits = rev) +
        scale_color_manual(expression(paste("total ", PM[2.5],"-trend classification")),
                           values = pmtrend_colors,
                           aesthetics = c("color", "fill")) +
        xlab("Number of states") + ylab("") + 
        scale_x_continuous(breaks = seq(0, 48, by = 8)) + 
        theme_classic() + 
        theme(legend.position = "bottom", 
              legend.direction = "vertical")}, 
  labels = "auto") %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureS07_pmtrend_sensitivity.pdf"), 
         ., width = 8, height = 7)

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

# table S3: smoke-influenced classes sensitivity ---- 
all_fits %>% 
  filter(state_code != "US") %>% 
  filter(PM_measure == "mean") %>%
  select(spec_factor, smoke_group, smoke_influenced_text) %>%
  mutate(smoke_group = gsub(",\n", ", ", smoke_group)) %>%
  {full_join(group_by(., spec_factor, smoke_group) %>% 
               summarise(n = n(), .groups = "drop") %>% 
               pivot_wider(names_from = smoke_group, values_from = n,
                           values_fill = 0), 
             group_by(., spec_factor, smoke_influenced_text) %>% 
               summarise(n = n(), .groups = "drop") %>% 
               pivot_wider(names_from = smoke_influenced_text, values_from = n,
                           values_fill = 0))} %>% 
  select(spec_factor, `no smoke influence detected`, `smoke-influenced`, 
         `smoke-caused reversal`, `smoke-influenced reversal`, 
         `smoke-influenced stagnation`, `smoke-influenced, no early decline`, 
         `smoke-influenced, no signficant trend change`) %>%
  arrange(spec_factor) %>% 
  mutate(spec_factor = as.character(spec_factor)) %>%
  stargazer::stargazer(type = "latex", summary = FALSE,
                       rownames = FALSE, 
                       font.size = "footnotesize",
                       # title = "\\textbf{Counts of states in different classifications under different sample restrictions and/or statistical specifications}. Main sample uses station-years that report at least 50 days in each year, with the year break in 2016. Other samples are as listed,  with \`\`yrs\" samples restricted to those reporting at least 50 days (unless otherwise specificed) in each that number of years; e.g. \`\`5yrs, 100obs\" restricts to stations that report at least 100 days in each of at least 5 years. \`\`Drop\" samples are those that drop individual years. \`\`No dup. break yr\" is a sample that does not duplicate the break year. \`\`Piecewise\" forces segments on either side of the break year to intersect at the break year.",
                       title = "\\textbf{Counts of states in different classifications under different sample restrictions and/or statistical specifications}. Samples and specifications are as in Table \\ref{table:pmtrends}.",
                       label = "table:smokeinfluence") %>% 
  magrittr::inset(11, "\\textit{Specification} & \\multicolumn{2}{c}{\\textbf{Smoke-influenced?}} & \\multicolumn{5}{c}{\\textbf{Influenced}} \\\\ \\cline{4-8}
\\textit{or Sample} & no & yes & caused & influenced  & influenced  & influenced, & influenced, \\\\ 
 &  &  & reversal & reversal &  stagnation & no early dec. & no sig. trend\\\\ ") %>% 
  writeLines(file.path(path_github, "tables", "raw", 
                       "tableS03_smokeinfluenced_sensitivity.tex"))

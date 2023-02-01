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

# extremes ----
extremePct_cdf <- station_year %>% 
  filter(totalPM_n > 50) %>% 
  group_by(id) %>% 
  mutate(n_year = n()) %>% 
  ungroup %>% 
  filter(n_year >= 15) %>% 
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
         group = ifelse(start == end, 
                        "annual", 
                        period)) 

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

# fig S12: extreme days from wildfire smoke by year/period ----
test <- station_year %>% 
  filter(totalPM_n > 50) %>% 
  group_by(id) %>% 
  mutate(n_year = n()) %>% 
  ungroup %>% 
  filter(n_year >= 15) %>% 
  mutate(pct_exceed = ifelse(totalPM_extremePct > 0, 
                             (totalPM_extremePct - nonsmokePM_extremePct)/
                               totalPM_extremePct, 
                             0)) %>% 
  {purrr::pmap_dfr(data.frame(start = c(2006, 2011, 2020, 2006:2022),
                              end = c(2010, 2022, 2022, 2006:2022)),
                   function(start, end){
                     filter(., year >= start & year <= end) %>%
                       # filter(., year >= 2006 & year <= 2010) %>%
                       group_by(year, state_abbr) %>% 
                       summarise(across(ends_with("exceed"), 
                                        mean),
                                 .groups = "drop") %>% 
                       group_by(state_abbr) %>% 
                       summarise(across(ends_with("exceed"), 
                                        mean),
                                 .groups = "drop") %>% 
                       arrange(desc(pct_exceed)) %>% 
                       # mutate(n_gte = 1:n()) %>% 
                       # group_by(pct_exceed) %>%
                       # summarise(n_gte = max(n_gte), 
                       #           states = list(state_abbr)) %>%
                       # # add row to get down to zero 
                       # rbind(., 
                       #       filter(., pct_exceed == max(pct_exceed)) %>% 
                       #         mutate(n_gte = 0)) %>%
                       mutate(start = start, 
                              end = end) %>% 
                       return
                   })} %>% 
  mutate(period = ifelse(start == end, 
                         start, 
                         paste0(start, "-", end)),
         type = ifelse(start == end, "annual", "range"),
         group = ifelse(start == end, period, "range")) 

plot_grid(
  # panel a
  test %>% 
    filter(type == "annual") %>%
    ggplot(aes(x = pct_exceed,
               y = period, 
               # height = ..ndensity..,
               group = period, 
               color = period,
               fill = period)) +
    geom_density_ridges(alpha = 0.6, color = alpha("grey20", 0.7),
                        # stat = "density", trim = TRUE, 
                        # scale = 1.4, 
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
  # panel b 
  extremePct_cdf %>%   
    filter(group != "annual" | period == "2021") %>% 
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
  ggsave(file.path(path_github, "figures", "raw", "figureS12_extremePct_from_smoke_v2.pdf"), 
         ., width = 7, height = 4.5)

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

# figs S3/S4: break point small multiples ----
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

purrr::pmap(data.frame(pm_measure = c("mean", "extremePct"), 
                       main_break = c(2016, 2010),
                       ylabel = I(c(expression(paste(PM[2.5]," (",mu, "g/", m^3, ")")), 
                                    expression(paste("% of days > 35 ",mu, "g/", m^3)))),
                       ylab_fn = I(c(NA, scales::percent)))[1,], 
            function(pm_measure, main_break, ylabel, ylab_fn){
              if (pm_measure == "mean") {
                fig_num = "S03"
              } else if (pm_measure == "extremePct") {
                fig_num = "S04"
              }
              
              if(is.na(ylab_fn)){
                y_lab_fn <- waiver()
              } else{
                y_lab_fn <- ylab_fn
              }
              
              fig <- breaks %>% 
                filter(PM_measure == pm_measure) %>% 
                left_join(state_dict) %>% 
                mutate(scale = "region") %>% 
                rbind(left_join(data.frame(scale = "CONUS", 
                                           PM_measure = pm_measure, 
                                           spec = "n50_y10",
                                           breakpoint = main_break, 
                                           nrow = NA), 
                                state_dict %>% 
                                  rbind(data.frame(state_abbr = "US", 
                                                   state_code = NA, 
                                                   state_name = "Contiguous US", 
                                                   climate_regions = NA)), 
                                by = character())) %>%
                filter(state_abbr != "DC") %>%
                {ggplot(data = .) + 
                    geom_rect(data = state_dict %>% 
                                filter(state_abbr != "DC"),
                              aes(fill = climate_regions),
                              xmin = -Inf,xmax = Inf,
                              ymin = -Inf,ymax = Inf, 
                              show.legend = F) +
                    geom_line(data = state_year %>% 
                                filter(pm_type == "totalPM") %>%
                                rename(value = all_of(pm_measure)),
                              aes(x = year, y = value,
                                  group = interaction(state_abbr, pm_type))) + 
                    # geom_vline(data = filter(., spec != "n50_y5"),
                    #            mapping = aes(xintercept = breakpoint), 
                    #            color = "grey30", lwd = 0.4) +
                    geom_vline(data = filter(., spec == "n50_y10"),
                               mapping = aes(xintercept = breakpoint, 
                                             linetype = scale), 
                               color = "black") +
                    facet_geo(~state_abbr, scales = "free_y",
                              grid = my_grid, label = "code") +
                    guides(linetype=guide_legend(title="Breakpoint scale")) +
                    theme_classic(base_size = 9.5) +
                    theme(strip.background = element_rect(linetype = "blank", fill = NA), 
                          strip.text = element_text(face = "bold", 
                                                    size = 10.5,
                                                    vjust = -0.5),
                          panel.spacing.x = unit(3, "points"),
                          panel.spacing.y = unit(-2, "points"),
                          legend.title = element_text(size = 13),
                          legend.position = c(0.89, 0.2),
                          legend.text = element_text(size = 12),
                          axis.title = element_text(size = 12.5),
                          plot.margin = unit(c(5.5, 8, 4, 5.5), "pt")) +
                    scale_color_manual(values = alpha(MetBrewer::met.brewer("Hiroshige", 9)[c(1, 7, 3:4, 2, 6, 5, 8:9)],
                                                      c(rep(1, 8), 0.8)),
                                       aesthetics = c("color", "fill")) + 
                    scale_x_continuous(breaks = c(2000, 2022),
                                       labels = c("   2000", "2022   ")) +
                    scale_y_continuous(breaks = scales::breaks_extended(n = 4, Q = c(1, 2, 3, 4, 5, 6)),
                                       labels = y_lab_fn) + 
                    xlab("") + ylab(ylabel)} 
              # make it a gtable
              fig %<>% ggplot_build %>% ggplot_gtable
              
              fig %<>% clean_geo_facet
              
              ggsave(file.path(path_github, "figures", "raw", 
                               paste0("figure", fig_num, "_breaks_", pm_measure, ".pdf")), 
                     fig, width = 11.5, height = 8.5)
            })

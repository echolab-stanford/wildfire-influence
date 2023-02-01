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

# fig 3/S10: mean and extreme pct small multiples ---- 
small_mult_figs <- purrr::pmap(
  data.frame(plot_spec = "main",
             pm_measure = c("mean", "extremePct"),
             grouping = c("smoke_group_factor", "smoke_influenced_text"),
             class_colors = I(list(alpha(c("#62205f", "#bb292c", "#ef8737", "#ffd353", "grey90"), 
                                         c(0.6, 0.7, 0.8, 0.8, 0.4)), 
                                   c(alpha("grey90", 0.4), "#98ab76"))), 
             ylabel = I(c(expression(paste(PM[2.5]," (",mu, "g/", m^3, ")")), 
                          expression(paste("% of days > 35 ",mu, "g/", m^3)))),
             ylab_fn = I(c(NA, scales::percent)), 
             fig_name = c("figure03_smoke_influence.pdf", 
                          "figureS10_extremePct_smoke_influence.pdf"))[1,], 
  function(plot_spec, pm_measure, grouping, class_colors, 
           ylabel, ylab_fn, fig_name){
    if(is.na(ylab_fn)){
      y_lab_fn <- waiver()
    } else{
      y_lab_fn <- ylab_fn
    }
    fig <- state_year %>% 
      rename(value = pm_measure) %>%
      ggplot() + 
      geom_rect(data = all_fits %>%
                  filter(spec == plot_spec & PM_measure == pm_measure) %>%
                  rename(color_by = all_of(grouping)),
                aes(fill = color_by),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf) +
      geom_line(aes(x = year, y = value,
                    group = interaction(state_abbr, pm_type), 
                    color = pm_type)) + 
      geom_vline(data = all_fits %>% 
                   filter(spec == plot_spec & PM_measure == pm_measure),
                 aes(xintercept = breakpoint), linetype = "dotted") +
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
      xlab("") + ylab(ylabel) 
    
    # make it a gtable
    fig %<>% ggplot_build %>% ggplot_gtable
    
    fig %<>% clean_geo_facet()
    
    ggsave(file.path(path_github, "figures", "raw", fig_name), 
           fig, width = 11.5, height = 8.5)
    
  })

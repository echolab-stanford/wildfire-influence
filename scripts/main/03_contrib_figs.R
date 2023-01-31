source("scripts/setup/01_load_packages.R")
source("scripts/setup/02_load_functions.R")
source("scripts/setup/03_load_settings.R")

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

# fig 2a: sankey of classifications ----
all_fits %>% 
  filter(PM_measure == "mean" & spec == "main" & state_abbr != "US") %>% 
  group_by(smoke_group, change_group) %>% 
  summarise(count = n(), .groups = "drop") %>%  
  gather_set_data(1:2) %>% 
  mutate(across(c(y, smoke_group, change_group),
                ~factor(.x,
                        levels = c("reversal",
                                   "stagnation",
                                   "no sig. early decline",
                                   "non-sig. change",
                                   "smoke-caused reversal",
                                   "smoke-influenced reversal",
                                   "smoke-influenced stagnation",
                                   "smoke-influenced,\nno early decline",
                                   "no smoke-influence detected"),
                        ordered = TRUE))) %>% 
  group_by(y) %>% 
  mutate(group_count = sum(count)) %>% 
  ungroup %>% 
  {ggplot(data = ., 
          aes(y = x, id = id, split = y, value = count)) +
      geom_parallel_sets(fill = "grey70",
                         #aes(fill = change_group), 
                         alpha = 0.3, axis.width = 0.1) +
      geom_parallel_sets_axes(aes(fill = y), 
                              axis.width = 0.1) +
      geom_text(data = mutate(., lab_position_y = ifelse(y == change_group, 2, 1)) %>% 
                  select(y, group_count, lab_position_y) %>% 
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
                  select(y, group_count, lab_position_y) %>% 
                  unique %>%
                  group_by(lab_position_y) %>% 
                  arrange(desc(y)) %>% 
                  mutate(lab_position_x = 2.5*((1:n()) - 1) + cumsum(group_count) - group_count/2) %>% 
                  mutate(axis_lab = gsub("^smoke-", "", y) %>% 
                           gsub(" stagnation", "\nstagnation", .) %>% 
                           gsub(" reversal", "\nreversal", .) %>% 
                           gsub("influence ", "influence\n",. )  %>% 
                           gsub("sig. early ", "early\n",. ) %>% 
                           gsub("sig.", "sig.\n", .)) %>% 
                  filter(lab_position_y == 1),
                aes(x = lab_position_x, 
                    y = lab_position_y, 
                    label = axis_lab), 
                vjust = 1, 
                size = 3,
                fontface = "italic",
                lineheight = 0.8,
                nudge_y = -0.08,
                inherit.aes = FALSE) +
      geom_text(data = mutate(., lab_position_y = ifelse(y == change_group, 2, 1)) %>% 
                  select(y, group_count, lab_position_y) %>% 
                  unique %>%
                  group_by(lab_position_y) %>% 
                  arrange(desc(y)) %>% 
                  mutate(lab_position_x = 2.5*((1:n()) - 1) + cumsum(group_count) - group_count/2) %>% 
                  mutate(axis_lab = gsub("^smoke-", "", y) %>% 
                           gsub(" stagnation", "\nstagnation", .) %>% 
                           gsub(" reversal", "\nreversal", .) %>% 
                           gsub("influence ", "influence\n",. )  %>% 
                           gsub("sig. early ", "early\n",. ) %>% 
                           gsub("sig.", "sig.\n", .)) %>% 
                  filter(lab_position_y == 2),
                aes(x = lab_position_x, 
                    y = lab_position_y, 
                    label = axis_lab), 
                vjust = 0, 
                fontface = "italic",
                nudge_y = 0.08,
                size = 3,
                inherit.aes = FALSE) +
      # geom_parallel_sets_labels(colour = 'black', angle = 0) + 
      # geom_parallel_sets_labels(aes(split = group_count), 
      #                           colour = 'black', angle = 0) + 
      scale_fill_manual(values = c(alpha("#175449", 0.85), # reversal
                                   "#6e948c", # stagnation 
                                   "#5d6174", # non-sig change
                                   "#dcdcdb", #alpha("#b9b9b8", 0.6), # no sig early decline 
                                   alpha("#62205f", 0.8), # smoke-caused revesal
                                   "#bb292c", # smoke-influenced reversal
                                   "#ef8737", # smoke-influenced stagnation
                                   "#ffd353", # smoke-influenced, no early decline
                                   "grey60")) + # not smoke influencned
      scale_y_continuous(expand = expansion(mult = 0.2)) +
      xlab("") + ylab("") + theme_void() + 
      theme(legend.position = "none")} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figure02_classification_sankey.pdf"), 
         ., width = 6, height = 2.5)

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
                          "figureS10_extremePct_smoke_influence.pdf")), 
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

# fig S1: smoke PM definition and station locations ----
# one panel with daily station data 
station_day <- readRDS(file.path(
  path_dropbox, "data", 
  "epa_station_day_totalPM_smokePM_panel_20000101-20221021_era5.rds"
)) %>% 
  filter(!is.na(pm25)) 
station_day %>% 
  filter(id == "060811001") %>%
  filter(date >= as.Date("2020-07-01") & date <= as.Date("2020-10-02")) %>%
  mutate(nonsmokePM = pm25 - smokePM, 
         totalPM = pm25) %>% 
  select(date, totalPM, nonsmokePM, smoke_day) %>% 
  mutate(nonsmoke_max = nonsmokePM, 
         nonsmoke_min = 0, 
         total_max = totalPM, 
         total_min = nonsmokePM) %>% 
  pivot_longer(ends_with(c('max', 'min')), 
               names_sep = "_", 
               names_to = c("pm_type", ".value")) %>%  
  mutate(pm_type = gsub("non", "non-", pm_type)) %>%
  arrange(pm_type) %>% 
  {ggplot(data = ., 
          aes(x = date)) + 
      geom_line(aes(y = max, color = pm_type, group = pm_type)) +
      geom_ribbon(aes(ymax = max, ymin = min, fill = pm_type), 
                  color = NA, alpha = 0.7) + 
      geom_point(aes(x = date, y = -3, 
                     alpha = I(ifelse(smoke_day, 1, 0))), 
                 color = "gray") + 
      scale_color_manual(values = c("blue", "black"),
                         aesthetics = c("fill", "color")) + 
      theme_classic() + 
      theme(legend.position = "none") + 
      xlab("") + 
      ylab(expression(paste("daily ", PM[2.5]," (",mu, "g/", m^3, ")")))} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureS01a_station-day_example.pdf"), 
         ., width = 3.5, height = 2.5)

# second with station-year non-smoke and smoke averages
station_year %>% 
  # filter(id == "060010009") %>%
  filter(id == "060811001") %>%
  pivot_longer(ends_with("mean")) %>% 
  mutate(name = gsub("PM_mean", "", name)) %>% 
  arrange(rev(name)) %>% 
  {ggplot(data = ., 
          aes(x = year, y = value, color = name, 
              linetype = name)) + 
      geom_line(lwd = 0.8) + 
      scale_color_manual(values = c("blue", "black"),
                         aesthetics = c("fill", "color")) + 
      scale_linetype_manual(values = c("32", "solid")) +
      theme_classic() + 
      theme(legend.position = "none") + 
      xlab("") + 
      ylab(expression(paste("annual ", PM[2.5]," (",mu, "g/", m^3, ")")))} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureS01b_station-year_example.pdf"), 
         ., width = 3.5, height = 2.5)

# lower panel with station locations by number of years with > 50 obs
station_ll <- read_sf(file.path(path_dropbox, "data", "epa_station_locations"))

states <- tigris::states(cb = TRUE, year = 2020) %>% 
  filter(STATEFP %in% nonContig_stateFIPS == FALSE)

station_year %>% 
  filter(totalPM_n > 50) %>% 
  group_by(id) %>% 
  summarise(n_year = n()) %>% 
  mutate(n_year = case_when(n_year < 5 ~ "< 5 years", 
                            n_year >= 5 & n_year < 10 ~ "5 - 9 years", 
                            n_year >= 10 & n_year < 15 ~ "10 - 14 years", 
                            n_year >= 15 ~ "15+ years"), 
         n_year = factor(n_year, 
                         levels = c("15+ years", 
                                    "10 - 14 years", 
                                    "5 - 9 years", 
                                    "< 5 years"), 
                         ordered = TRUE)) %>%
  left_join(station_ll %>% select(id = stn_id)) %>% 
  st_as_sf() %>% 
  st_transform(st_crs(states)) %>% 
  arrange(desc(n_year)) %>% 
  {ggplot(data = .) + 
      geom_sf(data = states, fill = NA) +
      geom_sf(aes(color = n_year), size = 0.8, alpha = 0.8) +
      scale_color_manual(name = "years with\nat least\n50 observations",
                         values = c("#591c19", "#b64f32", "#f7c267", "#8b8b99")) +
      theme_void() + 
      theme(legend.position = "right")} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureS01c_station_locations.pdf"), 
         ., width = 7, height = 4)

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
                       ylab_fn = I(c(NA, scales::percent))), 
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
# fig S5: coefficient estimate sensitivity ----
all_fits %>%
  filter(PM_measure == "mean") %>%
  select(state_name, t_df, spec, starts_with("b"), -breakpoint) %>% 
  pivot_longer(starts_with("b"), 
               names_sep = "_", 
               names_to = c("par", ".value")) %>% 
  filter(par %in% c("b1total", "b2total", "b2nonsmoke")) %>% 
  separate(par, into = c("period", "pm_type"), sep = 2) %>% 
  mutate(period = recode_factor(period, 
                                "b1" = "early period",
                                "b2" = "late period", 
                                .ordered = T), 
         pm_type = gsub("non", "non-", pm_type), 
         panel = paste0(period, ",\n", pm_type), 
         panel = factor(panel, 
                        levels = c("early period,\ntotal", 
                                   "late period,\ntotal", 
                                   "late period,\nnon-smoke"), 
                        ordered = TRUE),
         spec_type = case_when(spec == "main" ~ "main", 
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
      facet_wrap(~panel, scales = "free_x") + 
      scale_color_manual(values = c("#a82203", "#208cc0", "#f1af3a",
                                    "#cf5e4e", "#637b31", "black"),
                         name = "specification type",
                         guide = guide_legend(reverse = TRUE, 
                                              override.aes = list(shape = c(1, 16, 16, 16, 16, 16)))) +
      scale_y_discrete(limits = rev) +
      theme_classic() + 
      theme(strip.background = element_blank(), 
            strip.text = element_text(size = 11, face = "bold")) + 
      ylab("") + 
      xlab(expression(paste(PM[2.5]," trend (",mu, "g/", m^3, "/year)")))} %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureS05_coef_sensitivity.pdf"), 
         ., width = 7, height = 7)

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


# table S2: trend-change classes sensitivity ----
all_fits %>% 
  filter(state_code != "US") %>% 
  filter(PM_measure == "mean") %>% 
  group_by(spec_factor, change_group) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  pivot_wider(names_from = change_group, values_from = n,
              values_fill = 0) %>% 
  select(spec_factor, reversal, stagnation, acceleration, 
         `non-sig. change`, `no sig. early decline`) %>%
  arrange(spec_factor) %>% 
  mutate(spec_factor = as.character(spec_factor)) %>%
  stargazer::stargazer(type = "latex", summary = FALSE,
                       rownames = FALSE, 
                       font.size = "footnotesize",
                       title = "\\textbf{Counts of states in different total \\pmt-trend classifications under different sample restrictions and/or statistical specifications}. Main sample uses station-years that report at least 50 days in each year, with the year break in 2016. Other samples are as listed,  with \`\`yrs\" samples restricted to those reporting at least 50 days (unless otherwise specificed) in each that number of years; e.g. \`\`5yrs, 100obs\" restricts to stations that report at least 100 days in each of at least 5 years. \`\`Drop\" samples are those that drop individual years. \`\`No dup. break yr\" is a sample that does not duplicate the break year. \`\`Piecewise\" forces segments on either side of the break year to intersect at the break year.",
                       label = "table:pmtrends") %>% 
  magrittr::inset(11, "\\textit{Specification} \\\\ 
                  \\textit{or Sample} & reversal & stagnation & acceleration & non-sig. change & no sig. early decline \\\\ ") %>% 
  writeLines(file.path(path_github, "tables", "raw", 
                       "tableS02_pmtrend_sensitivity.tex"))

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

# fig S8: smoke-influenced classes sensitivity ---- 
smokeinfluence_colors <- c("#62205f", "#bb292c", "#ef8737", "#ffd353", "#dea868", "grey95")

plot_grid(
  # panel a
  all_fits %>% 
    filter(PM_measure == "mean") %>% 
    {ggplot(data = ., 
            aes(y = state_name,
                x = spec_factor, 
                color = smoke_group_factor, 
                fill = smoke_group_factor)) + 
        geom_tile() + 
        scale_y_discrete(limits = rev) + 
        scale_color_manual("Smoke-influence classification",
                           values = smokeinfluence_colors,
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
                color = smoke_group_factor, 
                fill = smoke_group_factor)) + 
        geom_bar() + 
        scale_y_discrete(limits = rev) +
        scale_color_manual("Smoke-influence classification",
                           values = smokeinfluence_colors,
                           aesthetics = c("color", "fill")) +
        xlab("Number of states") + ylab("") + 
        scale_x_continuous(breaks = seq(0, 48, by = 8)) + 
        theme_classic() + 
        theme(legend.position = "bottom", 
              legend.direction = "vertical")}, 
  labels = "auto") %>% 
  ggsave(file.path(path_github, "figures", "raw", "figureS08_smokeinfluenced_sensitivity.pdf"), 
         ., width = 8, height = 7)

# fig S9: western states sensitivity to break years -----
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

# table S4: extreme days sensitivity ----
all_fits %>% 
  filter(state_code != "US") %>% 
  filter(PM_measure == "extremePct") %>% 
  select(spec_factor, smoke_influenced_text) %>%
  group_by(spec_factor, smoke_influenced_text) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  pivot_wider(names_from = smoke_influenced_text, values_from = n,
              values_fill = 0) %>% 
  arrange(spec_factor) %>% 
  mutate(spec_factor = as.character(spec_factor)) %>%
  stargazer::stargazer(type = "latex", summary = FALSE,
                       rownames = FALSE, 
                       font.size = "footnotesize",
                       #title = "\\textbf{Counts of states in different extreme day smoke-influence classifications under different sample restrictions and/or statistical specifications}. Samples and specifications are as in Table \\ref{table:pmtrends}.",
                       label = "table:extremes_smokeinfluence") %>% 
  writeLines(file.path(path_github, "tables", "raw", 
                       "tableS04_extreme_smokeinfluenced_sensitivity.tex"))

# fig S11: extreme days sensitivity ----

plot_grid(
  # panel a
  all_fits %>% 
    filter(PM_measure == "extremePct") %>% 
    {ggplot(data = ., 
            aes(y = state_name,
                x = spec_factor, 
                color = smoke_influenced_text, 
                fill = smoke_influenced_text)) + 
        geom_tile() + 
        scale_y_discrete(limits = rev) + 
        scale_color_manual("Smoke-influence classification",
                           values = c("grey90", "#5a5a83"),
                           aesthetics = c("color", "fill")) +
        ylab("") + xlab("") + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1), 
              legend.position = "none")}, 
  # panel b
  all_fits %>% 
    filter(PM_measure == "extremePct") %>% 
    filter(state_abbr != "US") %>% 
    {ggplot(data = ., 
            aes(y = spec_factor, 
                color = smoke_influenced_text, 
                fill = smoke_influenced_text)) + 
        geom_bar() + 
        scale_y_discrete(limits = rev) +
        scale_color_manual("Smoke-influence classification",
                           values = c("grey90", "#5a5a83"),
                           aesthetics = c("color", "fill")) +
        xlab("Number of states") + ylab("") + 
        scale_x_continuous(breaks = seq(0, 48, by = 8)) + 
        theme_classic() + 
        theme(legend.position = "bottom", 
              legend.direction = "vertical")}, 
  labels = "auto") %>% 
  ggsave(file.path(path_github, "figures", "raw", 
                   "figureS11_extreme_smokeinfluenced_sensitivity.pdf"), 
         ., width = 8, height = 7)

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

# main text numbers ----
# abstract, 2/3 of states
all_fits %>% 
  filter(spec == "main" & PM_measure == "mean" & state_abbr != "US") %>% 
  group_by(smoke_influenced) %>% 
  summarise(pct = n()/48)

# abstract
all_fits %>% 
  filter(spec == "main" & PM_measure == "mean" & state_abbr != "US") %>% 
  filter(smoke_influenced) %>%
  transmute(state_abbr,
            b1total_est, 
            b2smoke_est = b2total_est - b2nonsmoke_est, 
            b1total_change = b1total_est*16, 
            b2smoke_change = (b2smoke_est)*6,
            pct_reversal = pmin(-b2smoke_change/b1total_change, 1),
            # pct_reversal = -b2smoke_change/b1total_change,
            years_undone = -b2smoke_change/b1total_est) %>% View
summarise(across(where(is.numeric), median), 
          pct_states = n()/48)

# number of extreme days in 2000 in eastern states 
station_year %>% 
  filter(year %in% 2000:2002) %>% 
  filter(climate_regions %in% c("Northeast", "Ohio Valley")) %>%
  summarise(avg_extremePct = mean(totalPM_extremePct)) %>% 
  multiply_by(365)

# impact in stagnating states
all_fits %>% 
  filter(spec == "main" & PM_measure == "mean" & state_abbr != "US") %>% 
  # filter(smoke_group == "smoke-influenced stagnation") %>%
  filter(smoke_influenced & !change_group == "no sig. early decline") %>% 
  transmute(state_abbr,
            change_group,
            b1total_est, 
            b2smoke_est = b2total_est - b2nonsmoke_est, 
            b1total_change = b1total_est*16, 
            b2smoke_change = (b2smoke_est)*6,
            pct_reversal = pmin(-b2smoke_change/b1total_change, 1),
            # pct_reversal = -b2smoke_change/b1total_change,
            years_undone = -b2smoke_change/b1total_est) %>% 
  group_by(change_group) %>% 
  summarise(across(where(is.numeric), median), 
            n = n())

all_fits %>% 
  filter(spec == "main" & PM_measure == "mean" & state_abbr != "US") %>% 
  # filter(smoke_group == "smoke-influenced stagnation") %>%
  filter(smoke_influenced) %>% 
  transmute(state_abbr,
            change_group,
            b1total_est, 
            b2smoke_est = b2total_est - b2nonsmoke_est, 
            b1total_change = b1total_est*16, 
            b2smoke_change = (b2smoke_est)*6,
            pct_reversal = pmin(-b2smoke_change/b1total_change, 1),
            # pct_reversal = -b2smoke_change/b1total_change,
            years_undone = -b2smoke_change/b1total_est) %>% 
  summarise(across(where(is.numeric), median), 
            n = n())

# extreme pct numnbers 
extremePct_cdf %>% 
  filter(pct_exceed >= 0.25) %>% 
  group_by(period) %>% 
  summarize(n = max(n_gte), 
            states = list(unique(states))) %>% 
  filter(grepl("-", period)) %>%
  View

extremePct_cdf %>% 
  filter(pct_exceed >= 0.75) %>% 
  filter(period == "2020-2022") %>%
  summarize(n = max(n_gte), 
            states = list(unique(states))) %>% 
  View

# methods numbers 
station_year$id %>% unique %>% length

feols(pm25 ~ smoke_day | id + state^month + state^year, 
      data = station_day %>% 
        mutate(year = lubridate::year(date),
               month = lubridate::month(date)))

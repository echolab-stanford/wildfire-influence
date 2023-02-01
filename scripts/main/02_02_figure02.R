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
                                   "no smoke influence detected"),
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
  ggsave(file.path(path_github, "figures", "raw", "figure02a_classification_sankey.pdf"), 
         ., width = 6, height = 2.5)

# fig 2b ----
#dt <- read_rds('data/state_pooled_fit_coefs_classifications.rds')
dt <- read_rds(file.path(path_dropbox, "output", "state_allFits_coefs_classifications.rds"))

dt <- dt %>%   mutate(b2smoke = b2total_est-b2nonsmoke_est, chg1=b1total_est*length(2001:2016),chg2=b2smoke*length(2017:2022), ratio=-chg2/chg1*100) %>%
  mutate(ratiot=ifelse(ratio>100,100,ratio))
df <- dt %>%  filter(spec=="main" & smoke_influenced==TRUE & state_abbr != "US" & PM_measure=="mean") 

pdf(file.path(path_github, "figures", "raw", "figure02b_PMchange.pdf"),width=5,height=3)
par(mfrow=c(2,1),mar=c(3,4,1,0))
hist(df$chg1,xlab="",las=1,xlim=c(-10,3),ylab="",main="",ylim = c(0,8),breaks=seq(-10,1,0.5),col="lightblue")
abline(v=median(df$chg1),col="lightblue",lty=2)

#hist(-df$chg2,xlim=xlim,las=1,ylab="",main="")
hist(df$chg2,add=T,main="",col="orange")
abline(v=median(df$chg2),col="orange",lty=2,lwd=2)

hist(df$ratiot,ylab="",xlab="",main="",las=1,breaks = seq(0,100,5))
abline(v=median(df$ratio),col="red",lty=2,lwd=2)
dev.off()

# ------------------------------------------------------------------------------
# Written by: Marshall Burke
# Plots Figure 4.
# ------------------------------------------------------------------------------
source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

# Load vpd data
historical_vpd <- read_rds(file.path(path_dropbox, 'data', 
                                     'gridmet_vpd_1979_2022_1deg_forest_cover.rds')) %>% 
  filter(west == 1)
vpd_obs_month <- historical_vpd %>% 
  filter(month > 4, month < 10) %>% 
  group_by(year, month) %>% 
  summarise(vpd = weighted.mean(vpd, forest, na.rm = T)) %>% 
  ungroup()
vpd_obs_year <- vpd_obs_month %>% 
  group_by(year) %>% 
  summarise(vpd = mean(vpd, na.rm = T)) %>% 
  ungroup() 

# Load smokePM
dt <- read_rds(file.path(path_dropbox, "data", "epa_station_year_2000_2022.rds"))
sc <- read_rds(file.path(path_dropbox, 'data', 'statecodes.rds'))
dt <- left_join(dt,sc)

# Define Western states
ws <- c("California", "Oregon", "Washington", "Nevada", "Arizona", "New Mexico", 
        "Utah", "Idaho", "Montana", "Wyoming", "Colorado")
dty <- dt %>% 
  mutate(smokePM = totalPM_mean-nonsmokePM_mean) %>% 
  filter(state %in% ws) %>% 
  group_by(year) %>% 
  summarise(smokePM = mean(smokePM,na.rm = T))

# Load gcm data
gcm <- read_rds(file.path(path_dropbox, 'data', 'CMIP6_west_vpd_debias_annual_2015_2100.rds'))
gcms <- gcm %>% 
  filter(year >= 2040 & year <= 2060) %>%
  group_by(model, scenario) %>% 
  summarise(vpd = mean(vpd_debias))

# Join data
toplot = full_join(dty, vpd_obs_year)

# ------------------------------------------------------------------------------
# Figure 4: climatic drivers of smoke ----
pdf(file = file.path(path_github, "figures", "raw", "figure04_VPD_smoke.pdf"), 
    height = 4, width = 8)

par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(toplot$vpd, log(toplot$smokePM), 
     pch = 21, axes = F, 
     xlab = "May-Sep avg VPD (kPa)", ylab = "annual smoke PM (ug/m3)", 
     bg = alpha("grey", 50), xlim = c(0.85, 1.1))
axis(1)
att = c(0.1, 0.2, 0.5, 1, 2, 3)
axis(2, las = 1, at = log(att), labels = att)
text(toplot$vpd + 0.02, log(toplot$smokePM), labels = toplot$year, cex = 0.7)
abline(lm(log(smokePM) ~ vpd,toplot), lty = 2)

plot(vpd_obs_year, type = "l", axes = F, xlim = c(1980, 2032), 
     ylim = c(0.75, 1.25), ylab = "May-Sep avg VPD (kPa)")
att = seq(1980, 2020, 10)
axis(1, at = att, att)
axis(2, las = 1)
abline(lm(vpd ~ year, vpd_obs_year), lty = 2)
ss <- c("ssp126", "ssp245", "ssp370")
colz = met.brewer("Hiroshige", 3, direction = -1)
xloc = c(2024, 2027, 2030)
w = 2  #width
for (s in 1:length(ss)) {
  out <- gcms %>% filter(scenario == ss[s]) 
  qq <- quantile(out$vpd, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  rect(xloc[s], qq[1], xloc[s] + w, qq[5], col = alpha(colz[s], 0.3), border = NA)
  rect(xloc[s], qq[2], xloc[s] + w, qq[4], col = alpha(colz[s], 0.6), border = NA)
  segments(xloc[s], qq[3], xloc[s] + w, qq[3], col = colz[s], lwd = 2)
}

dev.off()

# ------------------------------------------------------------------------------
# Sensitivity of vpd estimates to inclusion of time trend ----------------------
summary(lm(log(smokePM) ~ vpd, toplot))
summary(mod<- lm(log(smokePM) ~ vpd + year, toplot))


# Compute change in smoke at ensemble median VPD -------------------------------
ss <- c("ssp126", "ssp245", "ssp370")
meds <- gcms %>% 
  filter(scenario %in% ss) %>% 
  group_by(scenario) %>% 
  summarise(vpd = median(vpd))

# Add base period mean vpd to evaluate model
meds <- rbind(meds, 
              data.frame(scenario = "base", 
                         vpd = mean(toplot$vpd[toplot$year >= 2018])))

summary(mod <- lm(log(smokePM) ~ vpd, toplot))
s = summary(mod)$sigma
a0 = sum(exp(mod$residuals))/length(mod$residuals) # smearing estimate, eq 6.43 in wooldridge
meds <- meds %>% 
  mutate(pred = exp(predict(mod, meds))*exp(s^2/2), 
         predsmear = a0*exp(predict(mod, meds)))

# Two different ways of doing the predictions:
# First assumes normal distribution, 
# Second is smearing estimate that does not assume normality
# Second prediction more conservative in our case so we use that
data.frame(scen = meds$scenario[1:3], pred = meds$pred[1:3] - meds$pred[4])
data.frame(scen = meds$scenario[1:3], predsmear = meds$predsmear[1:3] - meds$predsmear[4])

# Distribution of smoke predictions over all models. 
# To get mean prediction, need to generate predictions and take mean rather than 
# compute mean VPD change and generate prediction because of nonlinearity
gcmp <- data.frame(gcms) %>% 
  mutate(pred = a0*exp(predict(mod, data.frame(vpd = gcms$vpd)))) %>% 
  filter(scenario != "ssp585") %>% 
  group_by(scenario) %>% 
  summarise(meanpred = mean(pred)) %>% 
  left_join(., meds)

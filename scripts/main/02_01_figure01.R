source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

f <- list.files(file.path(path_dropbox, "data", "epatrends"))
nms <- c("National","Northeast", "Northern Rockies","Northwest","Ohio Valley","South","Southeast","Southwest","Upper Midwest","West")
pdf(file=file.path(path_github, "figures", "raw", "figure01_PMTrends_smallmult.pdf"),width = 7,height=7)
par(mfrow=c(3,3),mar=c(3,3,1,1))
for (fl in 2:length(f)) {
  dt <- read_csv(file.path(path_dropbox, "data", "epatrends",f[fl]))
  plot(dt$Year,dt$Mean,axes=F,type = "l",xlab="",ylab="",col="grey")
  axis(1)
  axis(2,las=1)
  dt1 <- filter(dt,Year<=2016) %>% mutate(pred=predict(lm(Mean~Year)))
  lines(dt1$Year,dt1$pred,col="red")
  dt2 <- filter(dt,Year>=2016) %>% mutate(pred=predict(lm(Mean~Year)))
  lines(dt2$Year,dt2$pred,col="red")
  mtext(nms[fl],side=3,adj=1,cex=1,line=-2)
}
dev.off()

# national plot
pdf(file=file.path(path_github, "figures", "raw", "figure01_PMTrends_national.pdf"),width = 4,height=4)
par(mar=c(3,3,1,1))
dt <- read_csv(file.path(path_dropbox, "data", "epatrends", f[1]))
plot(dt$Year,dt$Mean,axes=F,type = "l",xlab="",ylab="",col="grey")
axis(1)
axis(2,las=1)
dt1 <- filter(dt,Year<=2016) %>% mutate(pred=predict(lm(Mean~Year)))
lines(dt1$Year,dt1$pred,col="red")
dt2 <- filter(dt,Year>=2016) %>% mutate(pred=predict(lm(Mean~Year)))
lines(dt2$Year,dt2$pred,col="red")
mtext("National",side=3,adj=1,cex=2,line=-2)
dev.off()


# look at trends robustness to outliers using Theil-Sen
dt <- read_csv(file.path(path_dropbox, "data", "epatrends",f[10]))
dt1 <- filter(dt,Year>2015) 
summary(lm(Mean ~ Year, dt1))
out <-  TheilSen(dt1$Year,dt1$Mean)

# write out state map
st <- states(cb=FALSE,resolution="20m") 
clim <- read_csv(file.path(path_dropbox, "data", "us_climate_regions.csv")) %>% rename(STUSPS=state_code)
st <- left_join(st,clim)
st <- st %>% filter(is.na(climate_regions)==F)
# pdf(file="output/statemap.pdf",width = 5,height=4)
jpeg(file=file.path(path_github, "figures", "raw", "figure01_statemap.jpg"),width = 1000,height=400)
gg <- ggplot()
gg <- gg + geom_sf(data=st,aes(fill=climate_regions),show.legend = F) + scale_fill_met_d("Hiroshige") + theme_minimal() + coord_sf(datum = NA)
gg
dev.off()

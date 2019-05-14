library(metR)
library(data.table)
library(ggplot2)
library(magrittr)


gh <- ReadNetCDF("~/DATOS/NCEP Reanalysis/hgt.daily.nc", 
                 c(gh = "hgt"),
                 subset = list(level = 300,
                               lat = -60:-45,
                               time = c("1979-01-01", "2018-12-31")))

gh <- gh[, .(gh = mean(gh)), by = .(time, lon)]

map <- geom_path(data = map_data("world2"), aes(long, lat, group = group))

ggplot() + 
  map +
  geom_vline(xintercept = c(217-90, 217+90))


filter_pacific <- function(lon, inverse = FALSE) {
  lon_range <- c(217-90, 217+90) # half hemisphere centered in 217
  
  ifelse(lon %between% lon_range, as.numeric(!inverse), as.numeric(inverse))
}

zw3_basins <- gh[, .(full_basin = FitWave(gh, 3)$amplitude,
                     pacific = 2*FitWave(gh*filter_pacific(lon, FALSE), 3)$amplitude,
                     no_pacific = 2*FitWave(gh*filter_pacific(lon, TRUE), 3)$amplitude),
                 by = .(time)]


zw3_basins %>% 
  # melt(id.vars = c("time", "full_basin"), variable.name = "basin") %>% 
  ggplot(aes(pacific, no_pacific)) +
  geom_point(size = 0.3, alpha = .1) +
  geom_smooth(method = "lm") +
  scale_y_continuous("Half hemisphere") +
  scale_x_continuous("Full hemisphere") +
  labs(title = "ZW3 in the Pacific is anticorrelated with ZW3 in the Atlantic + Indic")

zw3_basins[, ccf(pacific, no_pacific)]

t <- unique(gh$lon)
x <- cos(t*pi/180*3)

FitWave(x*filter_pacific(t, TRUE), 3)

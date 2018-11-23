library(metR)
library(data.table)
library(ggplot2)
library(magrittr)
library(ggperiodic)
library(here)
source(here("scripts/helperfun.R"))

map <- fortify(map_data("world2")) %>% 
   subset(lat < -20) %>% 
   geom_path(data = ., aes(long, lat, group = group), size = 0.3)

## Leo datos ----

hgt <- ReadNetCDF("~/DATOS/NCEP Reanalysis/hgt.mon.mean.nc", "hgt", 
           subset = list(level = 700, 
                         lat = c(-90, -20),
                         time = c("1979-01-01", "2000-12-31"))) %>% 
   setnames("level", "lev") %>% 
   .[, hgt.a := Anomaly(hgt), by = .(lon, lat, month(time))]

oficial <- fread("http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii") %>% 
   .[, time := lubridate::ymd_h(paste0(V1, "-", V2, "-01:00"))] %>% 
   .[, c("V1", "V2") := NULL] %>% 
   setnames("V3", "oficial")

# Dividiendo la AAO ----

copy(hgt) %>% 
   # .[, hgt := Anomaly(hgt), by = .(lon, lat, month(time))] %>%
   .[, aao := hgt.a*sqrt(cos(lat*pi/180))] %>%  
   EOF(aao ~ lon + lat | time, 1, data = .) -> aao

aao$left[, sim := mean(aao), by = .(lat)]
aao$left[, asim := Anomaly(aao), by = .(lat)]

ggplot(aao$left, aes(lon, lat)) +
   geom_contour_fill(aes(z = asim)) +
   geom_contour(aes(z = sim)) + 
   coord_polar()

hgt2 <- aao$left[, .(lon, lat, sim, asim)][hgt, on = c("lon", "lat")]


# Las aaos son los coeficientes de la regresión múltiple de la altura 
# geopotential y cada patrón. 

aaos <- ReadNetCDF("~/DATOS/NCEP Reanalysis/hgt.mon.mean.nc", "hgt", 
                  subset = list(level = 700, 
                                lat = c(-90, -20))) %>% 
   setnames("level", "lev") %>% 
   .[, hgt.a := Anomaly(hgt), by = .(lon, lat, month(time))] %>% 
   .[aao$left[, .(lon, lat, sim, asim)], on = c("lon", "lat")] %>% 
   .[, FitLm(-hgt.a, sim, asim), by = .(time)] %>% 
   .[term != "(Intercept)"] %>% 
   .[, estimate := (estimate - mean(estimate))/sd(estimate), by = term] 

oficial[aaos, on = "time"] %>% 
   ggplot(aes(oficial, estimate)) +
   geom_point(aes(color = term))

oficial[aaos, on = "time"] %>% 
   .[, cor(oficial, estimate, use = "complete.obs"), by = term]


# Efectos ----
aaos <- dcast(aaos, time ~ term, value.var = "estimate") 

temp <- ReadNetCDF(here::here("DATA/NCEP Reanalysis/air.mon.mean.nc"), c(t = "air"),
                   subset = list(level = 700, 
                                 lat = c(-90, -20),
                                 time = as.character(range(aaos$time)))) %>% 
   .[, t.a := Anomaly(t), by = .(lon, lat, month(time))]

aaos[temp, on = "time"] %>% 
   .[, FitLm(t.a, asim, sim), by = .(lon, lat)] %>%
   .[term != "(Intercept)"] %>% 
   ggperiodic::periodic(lon = c(0, 360)) %>% 
   # dcast(lon + lat ~ term, value.var = "estimate") %>% 
   ggplot(aes(lon, lat)) +
   geom_contour_fill(aes(z = estimate), breaks = AnchorBreaks(0, NULL, 0)) +
   map +
   scale_fill_divergent() +
   annotate("point", 310, -65) +
   coord_polar() +
   facet_wrap(~term)

(aaos[temp[lat %~% -65 & lon %~% 310], on = "time"] %>% 
   ggplot(aes(time, t.a)) +
   geom_line() +
   geom_line(aes(y = asim), color = "red") +
      geom_line(aes(y = sim), color = "blue")) %>% 
   
   ggwrap::ggwrap(2)

aaos[temp[lat %~% -65 & lon %~% 310], on = "time"] %>% 
   .[, cor := gtools::running(t.a, asim, cor, 11, pad = TRUE)] %>% 
   ggplot(aes(time, cor)) +
   geom_line() +
   geom_smooth()

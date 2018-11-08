library(metR)
library(data.table)
library(ggplot2)
library(magrittr)
library(ggperiodic)

hgt <- ReadNetCDF("~/DATOS/NCEP Reanalysis/hgt.mon.mean.nc", "hgt", 
           subset = list(level = 700, 
                         lat = c(-90, -20),
                         time = c("1979-01-01", "2000-12-31"))) %>% 
   setnames("level", "lev") %>% 
   .[, hgt.a := Anomaly(hgt), by = .(lon, lat, month(time))]

copy(hgt) %>%  
   .[, hgt := hgt.a*sqrt(cos(lat*pi/180)), by = .(lon, lat, month(time))] %>% 
   EOF(hgt ~ lon + lat | time, 1, data = .) -> aao
   
   
ggplot(aao$left, aes(lon, lat)) +
   geom_contour_fill(aes(z = hgt)) +
   coord_polar(start = 180*pi/180) +
   scale_fill_divergent() 

sam.pattern <- aao$right[, .(time, aao = hgt)] %>% 
   .[, aao := aao/sd(aao)] %>% 
   .[hgt, on = "time"] %>% 
   .[, FitLm(hgt, aao), by = .(lon, lat)] %>% 
   .[term == "aao"]

ggplot(sam.pattern, aes(lon, lat)) +
   geom_contour_fill(aes(z = -estimate), breaks = seq(-50, 20, by = 5)) +
   coord_polar(start = 180*pi/180) +
   scale_fill_divergent() 

oficial <- fread("http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii") %>% 
   .[, time := lubridate::ymd_h(paste0(V1, "-", V2, "-01:00"))] %>% 
   .[, c("V1", "V2") := NULL] %>% 
   setnames("V3", "oficial")




aao$right %>% 
   .[oficial, on = "time"] %>% 
   ggplot(aes(oficial, -hgt/sd(hgt, na.rm = TRUE))) +
   geom_point()

copy(hgt) %>% 
   # .[, hgt := Anomaly(hgt), by = .(lon, lat, month(time))] %>%
   .[, .(hgt = mean(hgt.a)*sqrt(cos(lat*pi/180))), by = .(lat, time)] %>%  
   EOF(hgt ~ lat | time, 1, data = .) -> zonal.eof

ggplot(zonal.eof$left, aes(lat, hgt)) +
   geom_line(aes(color = PC))

zonal.eof$right %>% 
   .[oficial, on = "time"] %>% 
   .[complete.cases(.)] %>% 
   ggplot(aes(-hgt/sd(hgt), oficial)) +
   geom_point()

zonal.pred <- predict(zonal.eof)

hgt[, mean(hgt.a), by = .(lat, time)] %>% 
   .[zonal.pred, on = c("lat", "time")] %>% 
   .[time == time[1]] %>% 
   melt(id.vars = c("lat", "time")) %>% 
   ggplot(aes(lat, value)) +
   geom_line(aes(color = variable))

zonal.pred[, .(lat, time, aaoz = hgt)] %>% 
   .[hgt, on = c("lat", "time")] %>%
   .[, hgt.aao := hgt.a - aaoz] -> hgt


copy(hgt) %>% 
   .[, hgt.aao := hgt.aao*sqrt(cos(lat*pi/180))] %>% 
   EOF(hgt.aao ~ lon + lat | time, 1:5, data = .) -> asim.eof

periodic(asim.eof$left, lon = c(0, 360)) %>% 
   ggplot(aes(lon, lat)) +
   geom_contour_fill(aes(z = -hgt.aao)) +
   facet_wrap(~PC) +
   coord_polar()

zonal.eof$right[, .(time, zonal.eof = hgt)] %>% 
   .[hgt, on = "time"] %>% 
   .[, cor(hgt.a, zonal.eof), by = .(lon, lat)] %>% 
   ggplot(aes(lon, lat)) +
   geom_contour_fill(aes(z = V1)) +
   scale_fill_divergent() +
   coord_polar()

sams <- sam.pattern[, .(lon, sam.zonal = mean(estimate),
                sam.asim = Anomaly(estimate)), by = .(lat)] 
ggplot(sams, aes(lon, lat)) +
   geom_contour_fill(aes(z = sam.asim)) +
   geom_contour2(aes(z = sam.zonal)) +
   coord_polar()

ggplot(sams, aes(lon, lat)) +
   geom_contour_fill(aes(z = sam.asim + sam.zonal)) +
   # geom_contour2(aes(z = sam.zonal)) +
   coord_polar()

sams.indexes <- sams[hgt, on = c("lat", "lon")] %>% 
   .[, .(zonal = cor(hgt, sam.zonal),
         asim = cor(hgt, sam.asim)), by = .(time)]

sams.indexes %>% 
   melt(id.vars = "time") %>% 
   ggplot(aes(time, value)) +
   geom_line(aes(color = variable))

sams.indexes %>% 
   .[oficial, on = "time"] %>% 
   .[complete.cases(.)] %>% 
   lm(oficial ~ zonal*asim, data = .) -> model



   melt(id.vars = "time") %>% 
   .[oficial, on = "time"] %>% 
   .[complete.cases(.)] %>% 
   
   .[, cor(value, oficial), by = variable]
   ggplot(aes(oficial, value)) +
   geom_point(aes(color = variable))
   

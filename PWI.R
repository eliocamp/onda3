library(metR)
library(data.table)
library(ggplot2)
library(lubridate)
library(magrittr)
library(here)

WaveEnvelope <- function(y, k = "all") {
   N <- length(y)
   x_hat <- fft(y)/N
   
   if (k[1] == "all") {
      k <- 1:ceiling(l/2)
   }
   
   x_hat[-k] <- 0
   Mod(fft(x_hat, inverse = T))*2
}

ncep <- ReadNetCDF(here("DATA/NCEP Reanalysis/vwnd.daily.nc"), c(v = "vwnd"),
                subset = list(time = c("1948-01-01", "2014-12-31"),
                              level = 500,
                              lat = -90:-20)) %>% 
   .[, gh := ReadNetCDF(here("DATA/NCEP Reanalysis/hgt.daily.nc"), c(gh = "hgt"),
                             subset = list(time = c("1948-01-01", "2014-12-31"),
                                           level = 500,
                                           lat = -90:-20), 
                        out = "vector")]

ncep[, v.s := RcppRoll::roll_mean(v, 31, fill = NA), by = .(lon, lat)]
ncep[!is.na(v.s), envelope := WaveEnvelope(v.s, 1:9), by = .(lat, time)]

PWI <- ncep[!is.na(v.s), max(envelope), by = .(lon, time)] %>% 
   .[, .(PWI = median(V1)), by = .(time)] 

copy(PWI) %>% 
   .[, PWI.s := smooth.loess(PWI ~ as.numeric(time), span = 365*5/.N)] %>%
   ggplot(aes(time, PWI)) +
   geom_line() +
   geom_line(aes(y = PWI.s), color = "red") 


wv <- copy(PWI) %>% 
   setnames("time", "date") %>% 
   WaveletComp::analyze.wavelet("PWI", make.pval = F)

abs(fft(PWI$PWI)) %>% 
   .[2:74] %>% 
   plot(type = "l")

sp <- PWI[, mean(PWI), by = .(yearmonth(time))] %>% 
   .[, PWI := Anomaly(V1), by = month(yearmonth)] %>% 
   .$PWI %>% 
   spectrum()

with(sp, data.table(freq, spec)) %>% 
   # .[1/freq < 400] %>% 
   ggplot(aes(1/freq, spec)) +
   geom_line() +
   scale_x_log10()

PWI %>% 
   ggplot(aes(month(time), PWI, group =month(time))) +
   geom_boxplot()
 
as.data.table(wv) %>% 
   ggplot(aes(date, period)) +
   geom_contour_fill(aes(z = amplitude)) +
   scale_y_continuous(trans = "log2")
 
ncep <- PWI[ncep, on = "time"]

ncep[!is.na(PWI), PWI.p := Percentile(PWI)]

means <- ncep[!is.na(PWI) & PWI >= .9, 
     .(gh = mean(gh)),
     by = .(lon, lat, season(time))] %>%
   .[, gh.z := Anomaly(gh), by = .(lat, season)] %>% 
   .[, ppsi := -gh.z*coriolis(-45)/coriolis(lat)]



ggplot(means, aes(lon, lat)) +
   geom_contour2(aes(z = ppsi, linetype = nsign(..level..))) +
   scale_y_latitude(limits = c(-90, -20)) +
   coord_polar() +
   facet_wrap(~season)


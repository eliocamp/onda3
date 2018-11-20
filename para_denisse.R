library(metR)
library(ggplot2)
library(data.table)
library(magrittr)

gh <- ReadNetCDF("DATA/NCEP Reanalysis/hgt.mon.mean.nc", c(gh = "hgt"),
                 subset = list(level = 200,
                              lat = -65:-40)) %>% 
   setnames("level", "lev") %>% 
   .[, .(gh = mean(gh)), by = .(lon, time)]

gh[, FitWave(gh, 3), by = .(time)] %>% 
   fwrite("lat_medias_mensual.txt")
   
gh <- ReadNetCDF("DATA/NCEP Reanalysis/hgt.daily.nc", c(gh = "hgt"),
                 subset = list(level = 200,
                               lat = -65:-40)) %>% 
   setnames("level", "lev") %>% 
   .[, .(gh = mean(gh)), by = .(lon, time)] %>% 
   .[, gh := RcppRoll::roll_mean(gh, 31, fill = NA)]

gh[!is.na(gh), FitWave(gh, 3), by = .(time)] %>% 
   fwrite("lat_medias_diaria.txt")

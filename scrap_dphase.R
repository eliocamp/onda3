library(metR)
library(data.table)
library(magrittr)
library(ggplot2)
library(circular)

gh <- ReadNetCDF("~/DATOS/NCEP Reanalysis/hgt.daily.nc", "hgt",
                 subset = list(lat = -55,
                               level = 300,
                               time = c("1980-01-01", "2018-01-01")))

qs <- gh[, FitWave(hgt, 1:4), by = .(lat, time)] %>% 
   .[, phase.c := circular(phase*k, modulo = "2pi")]



qs[, dphase := (phase - shift(phase))*k, by = k]
qs[dphase > pi, dphase := dphase - pi*2]
qs[dphase < -pi, dphase := dphase + pi*2]


qs %>% 
   .[, month := month(time)] %>% 
   .[!is.na(dphase)] %>% 
   ggplot(aes(phase*k, dphase)) +
   stat_summary_hex(aes(z = amplitude), 
                    bins = 10) +
   # geom_point(size = 0.1, alpha= 0.5) +
   geom_hline(data = qs[, weighted.mean(dphase, amplitude, na.rm = TRUE), 
                        by = k],
              aes(yintercept = V1, group = k)) +
   facet_wrap(~k) +
   scale_fill_viridis_c() 

qs %>% 
   ggplot(aes(abs(dphase), amplitude)) +
   geom_point(size = 0.3, alpha = 0.5) +
   scale_x_continuous(limits = c(0, pi))  +
   facet_wrap(~k)



library(metR)
library(data.table)
library(magrittr)


# Get SON mean geopotential height ----------------------------------------
# 
# 

hgt <- here::here("DATA", "hgt.mon.mean.nc") %>% 
  ReadNetCDF(vars = c("hgt"), 
             subset = list(
               lat = c(-90, -20),
               time = c("1981-01-01", "2010-12-31")), 
             key = TRUE) %>% 
  .[season(time) == "SON"] %>% 
  .[, .(hgt = mean(hgt)), by = .(lon, lat, time = seasonally(time))] %>% 
  .[, time := as.Date(time)]



# Compute complex eof of zonal anomalies ----------------------------------
eof <- hgt %>% 
  .[, hgt := Anomaly(hgt)*sqrt(cos(lat*pi/180)), by = .(time, lat)] %>% 
  .[, hgt := as.complex(hgt)] %>% 
  .[, hgt := spectral::analyticFunction(hgt), by = .(time, lat)] %>% 
  metR::EOF(hgt ~ time | lon + lat, data = ., n = 1:2) %>% 
  cut(2)
  

# Rotate eof to minimize Real correlation with ENSO -----------------------
rotate <- function(z, angle = 0) {
  complex(real = cos(angle), imaginary = sin(angle)) * z
}

enso <- rsoi::download_oni(TRUE, here::here("DATA", "oni.csv")) %>%
  as.data.table() %>%
  .[, .(time = as.Date(Date), oni = ONI)] %>%
  .[time %between% c("1981-01-01", "2010-12-31")] %>% 
  na.omit() %>%
  .[season(time) == "SON"] %>%
  .[, .(oni = mean(oni)), by = .(time = seasonally(time))]


with_enso <-  eof$left %>%
  copy() %>%
  .[enso, on = "time"] %>% 
  na.omit() 

angles <- seq(-pi, pi, by = .5*pi/180)

rotations <- lapply(angles, function(a) {
  with_enso %>%
    .[, hgt2 := rotate(hgt, a)] %>%
    .[, .(R = cor(Re(hgt2), oni),
          I = cor(Im(hgt2), oni))] %>%
    .[, rotation := a]
}) %>%
  rbindlist()

best_rotation <- rotations[I > 0][which.min(abs(R))]$rotation

eof$left[, hgt := rotate(hgt, best_rotation)]
eof$right[, hgt := rotate(hgt, best_rotation)]



# Extract spatial pattern -------------------------------------------------
pattern <- copy(eof$right)[, c("Real", "Imaginary") := list(Re(hgt), 
                                                         Im(hgt))] %>% 
  .[, hgt := NULL] %>% 
  .[, PC := NULL] %>% 
  melt(id.vars = c("lon", "lat"), variable.name = "part", value.name = "EOF")
  
  

# Write spatial pattern on disk -------------------------------------------
fwrite(pattern, here::here("DATA", "psa-pattern.csv"))


# Write timeseries --------------------------------------------------------
copy(eof$left)[, c("Real", "Imaginary") := list(Re(hgt), 
                                                 Im(hgt))] %>% 
  .[, hgt := NULL] %>% 
  .[, PC := NULL] %>% 
  melt(id.vars = c("time"), variable.name = "part", value.name = "EOF") %>% 
  fwrite(here::here("DATA", "psa-timerseries.csv"))

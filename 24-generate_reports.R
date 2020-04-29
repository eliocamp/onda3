library(magrittr)
library(data.table)
library(magrittr)
library(metR)

source("scripts/data_locations.R")
source("scripts/helperfun.R")

mean_season <- function(data) {
   data %>% 
      .[is.full_season(time)] %>%  
      .[, lapply(.SD, mean), by = .(lon, lat, time = seasonally(time))]   
}


with_cache <- function(file, fun) {
   file <- file.path("DATA", "cache", file)
   
   if (file.exists(file)) {
      return(readRDS(file))
   } else {
      out <- fun()
      saveRDS(out, file)
      return(out)
   }
}



subset <- list(lat = -90:0, 
               lev = list(50, 200))


files <- c(era   = ERA5(),
           # era20 = ERA20(),
           ncep  = NCEP())

var <- c("z", "hgt")
labs_datasets <- c(era20 = "ERA20C", era = "ERA5", ncep = "NCEP")
names(var) <- files

read_nc <- function(f) {
   nc <- ReadNetCDF(f, c(hgt = unname(var[f])), subset = subset) 
   
   if (unname(var[f]) == "z") {
      nc[, hgt := hgt/9.8]
   }
   return(nc)
}

hgt <- lapply(files, read_nc) %>% 
   rbindlist(idcol = "dataset") %>% 
   .[, mean_season(.SD), by = .(dataset, lev)] 

ks <- function(vorticity_gradient, u, lat) {
   vorticity_gradient/u*(metR:::a*cos(lat*pi/180))^2
}

sqrti <- function(x) {
   sqrt(abs(x))*sign(x)
}

datos <- with_cache("datos.Rds", function() {
   datos <- ReadNetCDF(ERA5(), vars = c(hgt = "z", "u", "v", vort = "vo", air = "t"),
                       subset = list(lat = c(-90:10))) %>% 
      na.omit() %>% 
      normalise_coords() %>% 
      .[, mean_season(.SD), by = lev]
   
   datos[, vort.dlat := Derivate(vort ~ lon + lat, cyclical = TRUE, sphere = TRUE)[[2]], 
         by = .(time, lev)] %>% 
      .[, U := mean(u), by = .(lat, lon, lev, season(time))] %>% 
      .[, ks := sqrti(ks(vort.dlat + f.dy(lat), U, lat))] %>% 
      .[, ":="(v_anom = Anomaly(v),
               u_anom = Anomaly(u)),
        by = .(lon, lat, lev, season(time))] %>%
      .[, ":="(v_star = Anomaly(v), 
               t_star = Anomaly(air)),
        by = .(lat, lev, time)] %>% 
      .[, ":="(vt_star = v_star*t_star,
               vu = v_anom*u_anom)] %>% 
      .[]
})


psi <- with_cache("psi.Rds", function() {
   ReadNetCDF(NCEP(vertical = "sigma"),
              subset = list(lat = c(-90, 80),
                            lev = c(0.2101),
                            time = c("1979-01-01", "2019-01-01"))) %>% 
      normalise_coords() %>%
      .[, mean_season(.SD)] %>% 
      
      .[, psi.z := Anomaly(psi), by = .(lat, lev, time)] %>% 
      .[, c("f.lon", "f.lat") := WaveFlux(.SD), by = .(lev, time)] %>% 
      .[]
}) 

file <- "20-capitulo_3.Rmd"
seasons <- c("SON", "DJF", "MAM", "JJA")
# seasons <- "SON"


purrr::walk(seasons, ~ rmarkdown::render(file, 
                                   params = list(season = .x, 
                                                 lats.eof = c(-90, -20)), 
                                   output_file = paste0("20b - ", .x)
                                   ))

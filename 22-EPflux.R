library(data.table)
library(magrittr)
library(ggplot2)
library(metR)

fecha <- "2002-10-01"
fecha <- rep(fecha, 2)

datos <- ReadNetCDF("~/DATOS/ERA-Interim/erai.mon.mean.nc", 
                    vars = c("v"),
                    subset = list(latitude = -90:90,
                                  time = fecha)) %>% 
  .[, u := ReadNetCDF("~/DATOS/ERA-Interim/erai.mon.mean.nc", 
                      vars = c("u"),
                      out = "vector",
                      subset = list(latitude = -90:90,
                                    time = fecha))] %>% 
  .[, t := ReadNetCDF("~/DATOS/ERA-Interim/erai.mon.mean.nc", 
                      vars = c("t"),
                      out = "vector",
                      subset = list(latitude = -90:90,
                                    time = fecha))]
a <-  6371000
datos %>% 
  .[, tita := Adiabat(level, t)] %>% 
  .[, logp := log(level*100)] %>% 
  .[,  dtita := 1/(level*100)*Derivate(tita ~ logp, equispaced = FALSE)[[1]], by = .(longitude, latitude, time)] 


datos %>% 
  .[, u_z := Anomaly(u), by = .(level, time, latitude)] %>% 
  .[, v_z := Anomaly(v), by = .(level, time, latitude)] %>% 
  .[, tita_z := Anomaly(tita), by = .(level, time, latitude)] %>% 
  .[, `:=`(v_zm = mean(v_z^2),
           u_zm = mean(u_z^2),
           uv = mean(u_z*v_z), 
           vtita = mean(v_z*tita_z)), by = .(level, time, latitude)] %>% 
  .[, `:=`(EPp = 1/dtita*coriolis(latitude)*a*tita_z * cos(latitude*pi/180),
           EPlon = 1/2*(v_zm - u_zm)*a*cos(latitude*pi/180)^2,
           EPlat = -a*uv*cos(latitude*pi/180)^2),
    by = .(longitude, latitude, time)] %>% 
  .[, ]


datos %>% 
  .[longitude %~% 180] %>% 
  .[time == time[1]] %>% 
  ggplot(aes(latitude, level)) +
  geom_arrow(aes(dx = EPlat*sqrt(1000/level)/(a*pi),
                 dy = EPp*sqrt(1000/level)*cos(latitude*pi/180)/1e5)) +
  # geom_streamline(aes(dx = Eplat, dy = EPp), L = 100) +
  # geom_arrow(aes(dx = Eplat,  dy = EPp)) +
  scale_mag() +
  coord_trans(y = metR::reverselog_trans())
scale_y_level()


datos %>% 
  # .[longitude %~% 180] %>% 
  .[level == 200] %>% 
  .[time == time[1]] %>% 
  .[is.cross(longitude, latitude)] %>% 
  ggplot(aes(longitude, latitude)) +
  geom_arrow(aes(dx = EPlon,
                 dy = EPlat)) +
  # geom_streamline(aes(dx = Eplat, dy = EPp), L = 100) +
  # geom_arrow(aes(dx = Eplat,  dy = EPp)) +
  scale_mag() 





fecha <- c(NA, NA)

datos <- ReadNetCDF("~/DATOS/ERA-Interim/erai.mon.mean.nc", 
                    vars = c("v"),
                    subset = list(latitude = 0:90,
                                  # level = 500,
                                  time = fecha)) %>% 
  .[, u := ReadNetCDF("~/DATOS/ERA-Interim/erai.mon.mean.nc", 
                      vars = c("u"),
                      out = "vector",
                      subset = list(latitude = 0:90,
                                    # level = 500,
                                    time = fecha))] %>% 
  .[, t := ReadNetCDF("~/DATOS/ERA-Interim/erai.mon.mean.nc", 
                      vars = c("t"),
                      out = "vector",
                      subset = list(latitude = 0:90,
                                    # level = 500,
                                    time = fecha))]

a <-  6371000
omega <- metR:::.omega
H <- 8000
datos %>% 
  .[, tita := Adiabat(level, t)] %>% 
  .[, dtp := metR:::.derv(t, level), by = .(longitude, latitude, time)] %>% 
  .[, `:=`(S = mean(-level*H*dtp + 2/7*t/H, na.rm = TRUE)), by = .(level, hemisphere = sign(latitude))] %>% 
  .[, `:=`(tita_z = Anomaly(tita),
           t_z = Anomaly(t),
           u_z = Anomaly(u),
           v_z = Anomaly(v)), 
    by = .(time, level, latitude)] %>% 
  .[, `:=`(vtita = v_z*tita_z,
           utita = u_z*tita_z,
           vt = v_z*t_z,
           uv = u_z*v_z,
           ttita = t_z*tita_z)] %>% 
  .[, `:=`(dvtita = metR:::.derv(vtita, longitude*pi/180, cyclical = TRUE),
           dutita = metR:::.derv(utita, longitude*pi/180, cyclical = TRUE),
           dttita = metR:::.derv(ttita, longitude*pi/180, cyclical = TRUE),
           p = level/1000)] %>% 
  # .[, `:=`(dtz = metR:::.derv(t, ))]
  .[, `:=`(Flat = p*cos(latitude*pi/180)*(v_z^2 - dvtita*2*omega*a*sin(2*latitude*pi/180)),
           Flon = p*cos(latitude*pi/180)*(-uv   + dutita*2*omega*a*sin(2*latitude*pi/180)),
           Fz   = p*cos(latitude*pi/180)/S*coriolis(latitude)*(vt  - dttita*2*omega*a*sin(2*latitude*pi/180)))] %>% 
  .[, `:=`(div = Divergence(Flon + Flat ~ longitude + latitude, cyclical = c(TRUE, FALSE))), 
    by = .(time, level)] %>% 
  .[, `:=`(divz = Divergence(Flat + Fz ~ latitude + level, cyclical = c(TRUE, FALSE))), 
    by = .(time, longitude)] 



datos %>% 
  .[level == 500] %>% 
  .[season(time) == "DJF"] %>% 
  .[, lapply(.SD, mean), by = .(longitude, latitude), .SDcols = c("Fz", "Flon", "Flat")] %>% 
  ggplot(aes(-longitude, -latitude)) +
  geom_contour(aes(z = Fz)) +
  geom_vector(aes(dx = Flon, dy = Flat), 
              data = function(d) d[is.cross(longitude, latitude)]) +
  scale_fill_divergent() +
  # geom_streamline(aes(dx = Flon, dy = Flat), L = 100) +
  scale_mag() +
  coord_polar()



datos %>% 
  .[season(time) == "DJF"] %>% 
  .[, lapply(.SD, mean), by = .(level, latitude), .SDcols = c("Fz", "Flon", "Flat")] %>% 
  ggplot(aes(latitude, level)) +
  # geom_contour_fill(aes(z = divz)) +
  geom_vector(aes(dx = Flat, dy = Fz*5e8)) +
  scale_fill_divergent() +
  scale_mag() +
  coord_trans(y = reverselog_trans())

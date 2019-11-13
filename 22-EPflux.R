library(data.table)
library(magrittr)
library(ggplot2)
library(metR)


fecha <- c("1998-11-01", "1998-11-01")

datos <- ReadNetCDF("~/DATOS/ERA-Interim/erai.mon.mean.nc", 
                    vars = c("v"),
                    subset = list(latitude = -90:90,
                                  # level = 500,
                                  time = fecha)) %>% 
  .[, u := ReadNetCDF("~/DATOS/ERA-Interim/erai.mon.mean.nc", 
                      vars = c("u"),
                      out = "vector",
                      subset = list(latitude = -90:90,
                                    # level = 500,
                                    time = fecha))] %>% 
  .[, t := ReadNetCDF("~/DATOS/ERA-Interim/erai.mon.mean.nc", 
                      vars = c("t"),
                      out = "vector",
                      subset = list(latitude = -90:90,
                                    # level = 500,
                                    time = fecha))]

datos %>% 
  .[, c("Flon", "Flat", "Flev") := EPflux(longitude, latitude, level, t, u, v), by = time]


datos %>% 
  .[level == 200] %>% 
  .[latitude <= 0] %>% 
  ggplot(aes(longitude, latitude)) +
  geom_contour_fill(aes(z = Flev)) +
  geom_vector(aes(dx = Flon, dy = Flat), 
              data =  function(d) d[is.cross(longitude, latitude, 0)]) +
  scale_fill_divergent() +
  # geom_streamline(aes(dx = Flon, dy = Flat), L = 100) +
  scale_mag() 




datos %>% 
  .[, .(Flat = mean(Flat), Flon = mean(Flon), Fz = mean(Fz)), 
    by = .(level, latitude)] %>%
  .[latitude <= 0] %>% 
  ggplot(aes(latitude, level)) +
  geom_contour_fill(aes(z = Fz), breaks = AnchorBreaks(0, NULL)) +
  geom_arrow(data = function(x) x[is.cross(latitude, level)],
             
             aes(dx = Flat*cos(latitude*pi/180)*sqrt(1000/level)/(a*pi),
                 dy = -Fz*cos(latitude*pi/180)*sqrt(1000/level)*2e2)) +
  geom_contour2(aes(z = Flat, linetype = factor(-sign(..level..)))) +
  # geom_vector(aes(dx = Flat, dy = Fz),
  # data = function(d) d[is.cross(longitude, latitude)]) +
  scale_fill_divergent() +
  scale_mag() +
  scale_y_level() +
  coord_cartesian(ylim = c(NA, 10))



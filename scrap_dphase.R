library(metR)
library(data.table)
library(magrittr)
library(ggplot2)
library(circular)
library(metR)

# from https://stackoverflow.com/questions/24198514/ggplot2-modify-geom-density2d-to-accept-weights-as-a-parameter
kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
   nx <- length(x)
   if (length(y) != nx) 
      stop("data vectors must be the same length")
   if (length(w) != nx & length(w) != 1)
      stop("weight vectors must be 1 or length of data")
   gx <- seq(lims[1], lims[2], length = n) # gridpoints x
   gy <- seq(lims[3], lims[4], length = n) # gridpoints y
   if (missing(h)) 
      h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y));
   if (missing(w)) 
      w <- numeric(nx)+1;
   h <- h/4
   ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
   ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
   z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
   # return(list(x = gx, y = gy, z = z))
   
   data.frame(expand.grid(x = gx, y = gy), z = as.vector(z))
}



gh <- ReadNetCDF("~/DATOS/NCEP Reanalysis/hgt.daily.nc", "hgt",
                 subset = list(lat = -55,
                               level = 300,
                               time = c("1980-01-01", "2018-01-01")))

qs <- gh[, FitWave(hgt, 1:4), by = .(time, lat)] %>% 
   .[, phase.c := circular(phase*k, modulo = "2pi")]


qs[, dphase := (phase - shift(phase))*k, by = k]
qs[dphase > pi, dphase := dphase - pi*2]
qs[dphase < -pi, dphase := dphase + pi*2]
qs[, dphase := dphase/k]


qs %>% 
   .[, month := month(time)] %>% 
   .[!is.na(dphase)] %>% 
   .[, kde2d.weighted(phase*k, dphase, amplitude), by = k] -> dens

qs %>% 
   .[, month := month(time)] %>% 
   .[!is.na(dphase)] %>%
   .[k == 3] %>% 
   ggplot(aes(dphase*k, phase*k)) +
   stat_summary_hex(aes(z = amplitude), 
                    bins = 20) +
   # geom_contour(data = dens, aes(x = x, y = y, z = z)) +
   # geom_point(size = 0.1, alpha= 0.5) +
   # geom_hline(data = qs[, weighted.mean(dphase, amplitude, na.rm = TRUE), 
   # by = k],
   # aes(yintercept = V1, group = k)) +
   facet_grid(season(month)~k) +
   scale_fill_viridis_c() 

qs %>% 
   .[, month := month(time)] %>% 
   .[!is.na(dphase)] %>% 
   # .[k == 3] %>% 
   ggplot(aes(abs(dphase)/pi, amplitude)) +
   geom_point(size = 0.3, alpha = 0.5) +
   # scale_x_continuous(limits = c(0, 360))  +
   facet_wrap(~k)


dates <- c("1980-01-01", "1980-02-28")
ks <- 2
gh %>% 
   .[, .(lon = lon, 
         hgt = FilterWave(hgt, ks)), by = time] %>% 
   .[time %between% dates] %>% 
   ggplot(aes(lon, time)) +
   geom_contour_fill(aes(z = hgt)) +
   geom_path(data = qs[time %between% dates & k == ks], 
             aes((phase[1] + cumsum(dphase))*180/pi, color = sign(dphase))) +
   scale_color_divergent()



qs[time %between% dates & k == ks] %>% 
   ggplot(aes(time, dphase)) +
   geom_line()

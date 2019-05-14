library(metR)
library(data.table)
library(ggplot2)
library(magrittr)


gh <- ReadNetCDF("~/DATOS/NCEP Reanalysis/hgt.daily.nc", 
                 c(gh = "hgt"),
                 subset = list(level = 300,
                               lat = -50:-30,
                               time = c("1979-01-01", "2018-12-31")))

gh <- gh[, .(gh = mean(gh)), by = .(time, lon)]

SH <- list(geom_path(data = subset(map_data("world2"), lat < 0),
                aes(long, lat, group = group), color = "gray50", size =0.3),
  geom_hline(yintercept = c(-50, -30), color = "gray50", size = 0.3))

ggplot_build.Layer <- function(plot) {
  plot <- ggplot() + 
    plot +
    theme_void()
  ggplot_build(plot)
}

ggplot_build.list <- function(plot) {
  plot <- ggplot() + 
    plot +
    theme_void()
  ggplot_build(plot)
}



PeriodicWavelet <- function(x, k) {
  period <- length(x)/k
  x1 <- rep(x, 3)
  keep <- (length(x)+1):(2*length(x))
  res <- list()
  for (p in seq_along(period)) {
    w <- WaveletComp::WaveletTransform(x1, dt = 1, upperPeriod = period[p],
                                       lowerPeriod = period[p])
    res[[paste0("k", ".", k[p])]] <- w$Ampl[keep]*sd(x)
    
  }
  return(res)
}

mblm.m <- memoise::memoise(mblm::mblm)
TheilSen <- function(formula, data, weights = NULL, prop = 1, repeated = TRUE) {
  sub <- withr::with_seed(42, sample(nrow(data), nrow(data)*prop, replace = FALSE))
  data <- data[sub, ]
  mblm.m(formula = formula, dataframe = data, repeated = repeated)
}



filter_pacific <- function(lon, inverse = FALSE) {
  lon_range <- c(217-90, 217+90) # half hemisphere centered in 217
  ifelse(lon %between% lon_range, as.numeric(!inverse), as.numeric(inverse))
}

zw3_basins <- gh[, .(full_basin = FitWave(gh, 3)$amplitude,
                     pacific = 2*FitWave(gh*filter_pacific(lon, FALSE), 3)$amplitude,
                     no_pacific = 2*FitWave(gh*filter_pacific(lon, TRUE), 3)$amplitude),
                 by = .(time)]

FitWave2 <- function(x, y, n) {
  sqrt(sum(FitLm(y, cos(x*n), sin(x*n))$estimate[-1]^2))
}

lon_range <- c(180, 360)
zw3_basins <- copy(gh) %>% 
  .[, gh := RcppRoll::roll_mean(gh, 90, fill = NA), by = .(lon)] %>%
  .[!is.na(gh), .(full_basin = FitWave(gh, 3)$amplitude,
                  western_hemisphere = FitWave2(lon[lon  %between% lon_range]*pi/180, gh[lon  %between% lon_range], 3),
                  eastern_hemisphere = FitWave2(lon[!(lon  %between% lon_range)]*pi/180, gh[!(lon  %between% lon_range)], 3)),
    by = .(time)]


zw3_basins %>% 
  .[year(time) == 1983] %>% 
  # melt(id.vars = c("time", "full_basin"), variable.name = "basin") %>% 
  ggplot(aes(western_hemisphere, eastern_hemisphere)) +
  stat_subset(geom = "point", aes(subset = time == min(time))) +
  geom_path() +
  geom_smooth(method = "lm") 

zw3_basins[, western_hemisphere - eastern_hemisphere] %>% 
  acf(lag.max = 30)

zw3_basins %>% 
  lm(full_basin ~ western_hemisphere + eastern_hemisphere, data  = .) %>% 
  summary()


zw3_basins %>% 
  melt(id.vars = c("time", "full_basin"), variable.name = "hemisphere") %>%
  ggplot(aes(full_basin, value)) +
  geom_point(size = 0.3, alpha = .1) +
  geom_smooth(method = "lm", formula = y ~ x - 1) +
  # scale_y_continuous("Half hemisphere") +
  # scale_x_continuous("Full hemisphere", trans = scales::log10_trans()) +
  facet_wrap(~reorder(hemisphere, -as.numeric(hemisphere)))
  # coord_trans(x = scales::exp_trans(10)) +
zw3_basins[, ccf(pacific, no_pacific, lag.max = 600)]


wavelet <- gh[TRUE][, gh := RcppRoll::roll_mean(gh, 31, fill = NA), by = .(lon)] %>% 
  .[!is.na(gh)] %>% 
  .[, .(lon = lon, amplitude = PeriodicWavelet(gh, 3)[[1]]), by = time]
  

wavelet[, FitWave(amplitude, 1), by = time] %>% 
  .[, density(circular::circular(phase, modulo = "2pi"), bw = 25)]  %>% 
  with(., data.table(x = x, y = y)) %>% 
  ggplot(aes(x*180/pi, y)) +
  annotation_custom(ggplotGrob(SH), ymax = .3, ymin = .1) +
  geom_line() +
  scale_y_continuous(limits = c(0, NA))
  
  
filter_wavelet <- function(lon, inverse = FALSE) {
  lon_range <- c(250-90, 250+90) # half hemisphere centered in 217
  ifelse(lon %between% lon_range, as.numeric(!inverse), as.numeric(inverse))
}

zw3_basins <- wavelet[, .(full_basin = mean(amplitude), 
                          pacific = mean(filter_wavelet(lon)*amplitude),
                          no_pacific = mean(filter_wavelet(lon, TRUE)*amplitude)), by = time]


wavelet[, mean(amplitude), by = lon] %>% 
  ggplot(aes(lon, V1)) + geom_line() +
  geom_vline(xintercept = c(250-45, 250+45, 70-45, 70+45))



wave_cor <- function(lon0) {
  wavelet[TRUE][, amplitude0 := amplitude[lon == lon0], by = time][, cor(amplitude0, amplitude), by = lon]
}

cors <- lapply(unique(wavelet$lon), wave_cor)
cors <- setNames(cors, unique(wavelet$lon))
wave_cors <- rbindlist(cors, idcol = "lon0")
wave_cors[, lon0 := as.numeric(lon0)]

ggplot(wave_cors, aes(lon, V1)) +
  geom_line() +
  gganimate::transition_manual(lon0)

mean_wavelet <- wavelet[, .(mean_amplitude = mean(amplitude)), by = time]

wavelet[, FitWave(amplitude, 1), by = time] %>%
  .[, dphase := (phase - shift(phase))*k, by = k] %>% 
  .[dphase > pi, dphase := dphase - pi*2] %>% 
  .[dphase < -pi, dphase := dphase + pi*2] %>% 
  .[mean_wavelet, on = "time"] %>% 
  ggplot(aes(amplitude/mean_amplitude, dphase*180/pi)) +
  geom_point(alpha = 0.1, size = 0.2) +
  scale_y_continuous("Desplazamiento (grados por día)") +
  scale_x_continuous("Amplitud") +
  facet_wrap(~k, scales = "free_y")

  
wavelet[, FitWave(amplitude, 1), by = time] %>%
  .[, acf(amplitude, lag.max = 300)]


wavelet[, c(FitWave(amplitude, 1),
            mean_amplitude = mean(amplitude)), by = time] %>% 
  .[, c("amplitude", "mean_mplitude") := .(Anomaly(amplitude), Anomaly(mean_amplitude)), by = yday(time)] %>% 
  .[, ccf(mean_amplitude, amplitude, lag.max = 300)]


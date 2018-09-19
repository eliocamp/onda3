



lats.index <- c(-65, -40)
levs.index <- c(100, 700)



ncep.day <- ReadNetCDF("DATA/hgt.daily_2000.nc", vars = c(gh = "hgt"),
                       subset = list(lat = lats.index, 
                                     level = levs.index))

index <- ncep.day[, qs3.index(gh, lat, level), by = time]



index %>% copy() %>% 
   .[, phase := phase*3] %>% 
   ggplot(aes(phase)) +
   
   geom_raster(stat = "density", aes(x = phase, y = yday(time), fill = ..density..,
                                     group = yday(time)))

+
   map.arg + 
   map.nz +
   scale_y_continuous("Month", breaks = 1:12, labels = month.abb, trans = "reverse",
                      minor_breaks = NULL) +
   scale_x_continuous("Phase", labels = LabDegrees, breaks = as.radians(seq(0, 360, by = 30)),
                      minor_breaks = NULL) +
   scale_fill_viridis_c(guide = "none") +
   coord_cartesian(ylim = c(1:12))



index[, .(amoma = stationarity.wave2(as.wave(amplitude, phase, 3))), 
      by = yday(time)] %>% 
   ggperiodic::periodic(yday = c(1, 367)) -> day_stat

ggplot(day_stat, aes(yday, amoma)) +
   geom_line() +
   geom_smooth(span = 0.2, n = 367)

acf(day_stat$amoma)




index[, .(amoma = stationarity.wave2(as.wave(amplitude, phase, 3))), 
      by = .(yearmonth(time))] -> monthly_stat

ggplot(monthly_stat, aes(yearmonth, amoma)) +
   geom_line() +
   geom_smooth(span = 0.2, n = nrow(monthly_stat))

acf.sig(monthly_stat$amoma, 30) %>% 
   no(1) %>% 
   ggplot(aes(lag, acf)) +
   geom_col(aes(fill = factor(sign(acf)))) 


index %>% 
   copy() %>% 
   .[, wave := as.wave(amplitude, phase, rep(3, .N))] -> indexx

index[, .(amoma = stationarity.wave2(as.wave(amplitude, phase, rep(3, .N)))), 
      by = .(isoweek(time))]  -> weekly

ggplot(weekly, aes(isoweek, amoma)) +
   geom_line() 

acf.sig(monthly_stat$amoma, 30) %>% 
   no(1) %>% 
   ggplot(aes(lag, acf)) +
   geom_col(aes(fill = factor(sign(acf)))) 

width <- 15

index[, amoma := listapply(as.wave(amplitude, phase, rep(3, .N)), width, stationarity.wave)]

acf.sig(index[!is.na(amoma), amoma], 30) %>% 
   no(1) %>% 
   ggplot(aes(lag, acf)) +
   geom_point()



ncep.day <- ReadNetCDF("DATA/hgt.daily_2000.nc", vars = c(gh = "hgt"),
                       subset = list(level = 200))
index[, .(amoma, time)][ncep.day, on = "time"] -> ncep.day

ncep.day[level == 200 & is.finite(amoma), FitLm(gh, amoma), 
         by = .(lon, lat, month(time))] %>% 
   .[term != "(Intercept)"] %>% 
   ggplot(aes(lon, lat)) +
   geom_contour_fill(aes(z = estimate), breaks = AnchorBreaks(0, NULL, 0)) +
   map.SH +
   scale_fill_divergent() +
   facet_wrap(~month)



qs2.index <- function(gh, lat, lev, lats.index =  c(-65, -40), levs.index = c(100, 700)) {
   dt <- data.table(gh, lat, lev)
   dt[lat %between% lats.index & 
         lev %between% levs.index] %>% 
      .[, FitWave(gh, 2), by = .(lat, lev)] %>% 
      .[, phase := circular(phase*2, modulo = "2pi")] %>% 
      .[, .(amplitude = mean(amplitude), phase = as.numeric(mean.circular(phase)/2))]
}

index2 <- ncep.day[, qs2.index(gh, lat, level), by = time]

index2 %>% copy() %>% 
   .[, phase := phase*2] %>% 
   ggplot(aes(phase)) +
   geom_raster(stat = "density", aes(x = phase, y = month(time), fill = ..density..,
                                     group = month(time)))


index2[, amoma := listapply(as.wave(amplitude, phase, rep(2, .N)), width, stationarity.wave)]
ncep.day[, amoma := NULL]
index2[, .(amoma, time)][ncep.day, on = "time"] -> ncep.day

ncep.day[level == 200 & is.finite(amoma), FitLm(gh, amoma), 
         by = .(lon, lat, month(time))] %>% 
   .[term != "(Intercept)"] %>% 
   ggplot(aes(lon, lat)) +
   geom_contour_fill(aes(z = estimate), breaks = AnchorBreaks(0, NULL, 0)) +
   map.SH +
   scale_fill_divergent() +
   facet_wrap(~month)

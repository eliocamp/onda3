smooths <- c(5, 10, 15, 30, 60, 90, 180, 365, 2*365 + 1, 3*365, 5*365, 
             10*365 +1 , 20*365 + 1 , 30*365 + 1)

v_spect <- function(data, col, s, ks = 1:9) {
   data[][, col.s := RcppRoll::roll_mean(get(col), s, fill = NA), by = .(lon)] %>% 
      .[!is.na(col.s)] %>% 
      .[, FitWave(col.s, ks), by = time] %>% 
      .[, phase := circular::circular(phase*k)] %>% 
      .[, lapply(.SD, mean), by = .(k), .SDcols = -"time"] 
}


ReadNetCDF("DATA/NCEP Reanalysis/hgt.daily.nc", c(gh = "hgt"), 
           subset = list(level = 500, 
                         lat = -55,
                         time = c("1979-01-01", "2014-12-31"))) -> gh

gh[, gh.a := Anomaly(gh), by = .(lon, lat, yday(time))]

data.table(s = smooths) %>% 
   .[, v_spect(gh, "gh", s, 1:4), by = s] -> perds_gh

data.table(s = smooths) %>% 
   .[, v_spect(gh, "gh.a", s, 1:4), by = s] -> perds_gh.a



ggplot(perds_gh, aes(s/31, r2)) +
   geom_line(aes(color = factor(k))) + 
   geom_point() +
   directlabels::geom_dl(aes(label = k), method = "last.points") +
   scale_x_log10("Smoothing window (months)", 
                 sec.axis = sec_axis(~.*31/365, "Smoothing window (years)")) +
   # scale_y_log10() +
   scale_color_viridis_d()


ggplot(perds_gh.a, aes(s/31, phase)) +
   geom_line(aes(color = factor(k))) + 
   geom_point() +
   directlabels::geom_dl(aes(label = k), method = "first.points") +
   scale_x_log10("Smoothing window (months)", 
                 sec.axis = sec_axis(~.*31/365, "Smoothing window (years)")) +
   # scale_y_log10() +
   scale_color_viridis_d()


ggplot(perds_gh, aes(s/31, amplitude, group = k)) +
   geom_line() +
   geom_point() +
   directlabels::geom_dl(aes(label = k), method = "last.points") +
   scale_x_log10("Smoothing window (months)", 
                 sec.axis = sec_axis(~.*31/365, "Smoothing window (years)")) +
   scale_y_log10() +
   scale_color_viridis_d() 

ggplot(perds_gh, aes(s/31, phase)) +
   geom_line(aes(color = factor(k))) +
   directlabels::geom_dl(aes(label = k), method = "last.points") +
   scale_x_log10("Smoothing window (months)", 
                 sec.axis = sec_axis(~.*31/365, "Smoothing window (years)")) +
   scale_color_viridis_d()  


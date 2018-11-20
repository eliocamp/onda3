smooths <- c(5, 10, 15, 30, 60, 90, 180, 365, 2*365 + 1, 3*365, 5*365, 
             10*365 +1 , 20*365 + 1 , 30*365 + 1)

v_spect <- function(data, col, s, ks = 1:9) {
   data[][, col.s := RcppRoll::roll_mean(get(col), s, fill = NA), by = .(lon)] %>% 
      .[!is.na(col.s)] %>% 
      .[, FitWave(col.s, ks), by = time] %>% 
      .[, phase := circular::circular(phase*k)] %>% 
      .[, lapply(.SD, mean), by = .(k, month(time)), .SDcols = -"time"] 
}


ReadNetCDF("DATA/NCEP Reanalysis/hgt.daily.nc", c(gh = "hgt"), 
           subset = list(level = 500, 
                         lat = -55,
                         time = c("1979-01-01", "2014-12-31"))) -> gh


data.table(s = smooths) %>% 
   .[, v_spect(gh, "gh", s, 1:20), by = s] -> perds_gh


ggplot(perds_gh, aes(s/31, r2)) +
   geom_line(aes(color = factor(k))) + 
   directlabels::geom_dl(aes(label = k), method = "last.points") +
   scale_x_log10() +
   scale_y_log10() +
   scale_color_viridis_d() +
   facet_wrap(~month)


ggplot(perds_gh, aes(s/31, amplitude)) +
   geom_line(aes(color = factor(k))) +
   directlabels::geom_dl(aes(label = k), method = "last.points") +
   scale_x_log10() +
   scale_y_log10() +
   scale_color_viridis_d() +
   facet_wrap(~month)

ggplot(perds_gh[k <= 5], aes(s/31, phase)) +
   geom_line(aes(color = factor(k))) +
   directlabels::geom_dl(aes(label = k), method = "last.points") +
   scale_x_log10() +
   scale_color_viridis_d()  +
   facet_wrap(~month)


perds_gh %>% 
   .[, sd(phase), by = .(month, k)] %>% 
   ggplot(aes(k, V1)) +
   geom_line() +
   facet_wrap(~month)


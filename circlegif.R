
lags <- -60:60
qs.eof$left %>% 
   dcast(time ~ PC, value.var = "gh") %>% 
   .[qs.trim(time) == "JFM"] %>% 
   {
      PC1 <<- .$PC1
      PC2 <<- .$PC2
   }
gh1 <- gh[time %in% unique(qs.eof$left$time)] %>% 
   .[qs.trim(time) == "JFM"] %>% 
   .[, gh.a := Anomaly(gh), by = .(lat, lon, yday(time))]



prediction <- function(x, x1, x2) {
   co <- x$coefficients
   new <- cbind(1, x1, x2)
   apply(new, 1, function(x) sum(x*co))
}



N <- 60
time <- 15  # total gif time in secons
bins <- 20

theta <- seq(0, 360, length.out = N)
if (theta[N] == 360) theta <- theta[-N]
M <- Mag(mean(PC1), mean(PC2))
x1 <- cos(theta*pi/180)*M
x2 <- sin(theta*pi/180)*M

gh1[, .(theta = theta, 
        gh.pred = prediction(.lm.fit(cbind(1, PC1, PC2),
                                     gh.a),
                             x1, x2)),
    keyby = .(lon, lat)] -> pred


binwidth <- pretty(diff(range(pred$gh.pred))/bins)[1]  # use pretty binwidth
ggplot(pred, aes(lon, lat)) +
   geom_contour_fill(aes(z = gh.pred, frame = factor(theta)), 
                     circular = "x", 
                     breaks = AnchorBreaks(0, binwidth, 0)) +
   map.SH +
   geom_vline(aes(xintercept = theta, frame = factor(theta))) +
   scale_fill_divergent(guide = guide_colorstrip_bottom(),
                        breaks = AnchorBreaks(0, binwidth, 0)) +
   scale_s_map(ylim = c(-90, -30)) +
   coord_polar() -> g 
# coord_quickmap() -> g

interval <- time/N
gganimate::gganimate(g, "lags.gif", interval = interval, title_frame = FALSE,
                     ani.width = 600, ani.height = 600)   



pred[theta %~% c(45, 135, 225, 315)] %>% 
   ggplot(aes(lon, lat)) +
   geom_contour_fill(aes(z = gh.pred, frame = factor(theta)), 
                     circular = "x", 
                     breaks = AnchorBreaks(0, binwidth, 0)) +
   map.SH +
   geom_vline(aes(xintercept = theta, frame = factor(theta))) +
   scale_fill_divergent(guide = guide_colorstrip_bottom(),
                        breaks = AnchorBreaks(0, binwidth, 0)) +
   scale_s_map(ylim = c(-90, -30)) +
   coord_polar() +
   facet_wrap(~theta)
   

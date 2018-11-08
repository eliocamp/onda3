library(data.table)

x <- rnorm(12*(2018-1948 +1))

r.s <- RcppRoll::roll_mean(x, 12*10, fill = NA)

wv <- Wavelets(data.frame(x = r.s[!is.na(r.s)]), "x")
wv <- as.data.table(wv)

wv[, upper := wavelet.boot(.signal.random, B = 1000)]

ggplot(wv, aes(location, period)) +
   geom_raster(aes(fill = amplitude)) +
   stat_subset(aes(subset = amplitude > upper), geom = "point", 
               size = 0.1, alpha = 0.4) +
   stat_subset(aes(subset = p.value < 0.05), geom = "point",
               size = 0.1, alpha = 0.4, color = "red") +
   scale_y_continuous(trans = "log2")

library(ggplot2)
library(metR)
library(data.table)
library(magrittr)

world <- subset(map_data("world2"), lat > 0)

geom_world <- geom_path(data = world, aes(long, lat, group = group),
                        size = 0.2)

data <- ReadNetCDF("~/DATOS/NCEP Reanalysis/uwnd.mon.mean.nc", c(var = "uwnd"),
                   subset = list(level = 1000,
                                 lat = c(-90:90),
                                 time = c("1948-01-01", "2009-12-31"))) %>% 
   .[month(time) %in% c(12, 1, 2)] %>% 
   .[, .(var = mean(var)), by = .(lat, lon, year(time))] 

random <- rnorm(uniqueN(data$year))

point <- data[data[lon %~% 100 & lat %~% 25, .(year, var_point = var)],
              on = "year"] %>% 
   .[, FitLm(var, var_point = year, se = TRUE), by = .(lon, lat)] %>% 
   .[, p.val := pt(-abs(estimate)/std.error, df)*2]


data[, linear := var - .lm.fit(cbind(1, year), var)$residuals, by = .(lon, lat)]
data[, constant := sign(linear[.N]) == sign(linear[1]), by = .(lon, lat)]

FDR <- function(p.vals, alpha = 0.01) {
   M <- length(p.vals)
   id <- order(p.vals)
   p.vals_ordered <- p.vals[id]
   
   a <- alpha * seq_along(p.vals)/M
   signif <- p.vals_ordered <= a
   signif[order(id)]
}


point[term == "var_point"] %>%
   .[lat %between% c(0, 90)] %>%
   # .[, signif := FDR(p.val)] %>%  
   ggplot(aes(lon, lat)) +
   geom_raster(aes(fill = estimate)) +
   # geom_contour_fill(aes(z = estimate)) +
   stat_subset(aes(subset = p.val < 0.01),
               geom = "point", size = 0.1, alpha = 0.3) +
   # stat_subset(aes(subset = p.adjust(p.val, "BY") <= 0.01),
               # geom = "point", size = 0.1, alpha = 0.3) +
   geom_world +
   scale_fill_divergent("Trend in U(300hpa)") +
   coord_quickmap() +
   metR:::theme_field() +
   annotate("point", x = 120, y = 70)


ggplot(data, aes(lon, lat)) +
   geom_raster(aes(fill = constant)) +
   # geom_contour_fill(aes(z = estimate)) +
   # stat_subset(aes(subset = p.val < 0.01),
               # geom = "point", size = 0.1, alpha = 0.3) +
   # stat_subset(aes(subset = p.adjust(p.val, "BY") <= 0.01),
   # geom = "point", size = 0.1, alpha = 0.3) +
   geom_world +
   # scale_fill_divergent("Trend in U(300hpa)") +
   coord_quickmap() +
   metR:::theme_field() +
   annotate("point", x = 120, y = 70)


ggplot(point, aes(p.val, p.adjust(p.val, "BY"))) +
   geom_point()

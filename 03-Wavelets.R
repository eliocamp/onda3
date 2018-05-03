---
   title: "02 - figuras para un paper"
author: "Elio Campitelli"
output:
   ioslides_presentation:
   fig_height: 5.5
fig_width: 10.5
smaller: yes
widescreen: yes
cache: yes
---
   
   ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE,
                      cache = TRUE, cache.lazy = TRUE,
                      warning = FALSE)

library(data.table)
library(ggplot2)
library(dplyr)
library(metR) 
library(WaveletComp)
# library(patchwork)
library(circular)
library(hrbrthemes)
library(extrafont)

source("scripts/helperfun.R")
# Plot thingys

data.world <- BuildMap(res = 1, smooth = 1)
map.world <- geom_map2(data.world)
map.SH <- geom_map2(data.world[lat %b% c(-90, 20)], color = "gray20")

# options(OutDec = ",")

pres <- ReadNetCDF("DATA/srfp.mon.mean.nc")
pres.mean <- pres[, .(pres = median(pres)), by = lat]
pres.mean <- rbind(data.table(lat = 0.0, pres = Inf), 
                   pres.mean, 
                   data.table(lat = -90.0, pres = Inf))
surface <- geom_polygon(data = pres.mean, aes(y = pres), fill = "white", 
                        alpha = 1, color = "gray30", size = 0.5)
pres <- pres[, .(pres = mean(pres)), by = .(lon, lat)]

# From https://github.com/hrbrmstr/hrbrthemes/issues/18
d <- read.csv(extrafont:::fonttable_file(), stringsAsFactors = FALSE)
d[grepl("Light", d$FontName),]$FamilyName <- font_rc_light
write.csv(d,extrafont:::fonttable_file(), row.names = FALSE)
extrafont::loadfonts()

theme_elio <- theme_minimal(base_size = 11) +
   theme(
      text = element_text(family = font_rc),
      legend.position = "bottom", legend.box = "vertical",
      panel.spacing.y = unit(5, "mm"),
      panel.spacing.x = unit(5, "mm"),
      legend.spacing = unit(2, "mm"),
      plot.margin = grid::unit(rep(3, 4), "mm"),
      legend.title = element_blank(),
      legend.box.spacing = unit(3, "mm"),
      legend.margin = margin(t = -5),
      panel.grid = element_line(color = "gray50", size = 0.2, linetype = 3),
      panel.ontop = TRUE)
theme_set(theme_elio)
update_geom_defaults(metR:::GeomTextContour, list(family = font_rc))

update_geom_defaults("contour2", list(color = "black"))
update_stat_defaults("contour2", aes(linetype = factor(-sign(..level..))))

coord_quickmap <- function(..., ylim = c(-90, 20)) {
   ggplot2::coord_quickmap(ylim = ylim, ...) 
} 

geom_contour_ <- function(..., gap = 0, rotate = FALSE) {
   if (gap != 0) {
      list(geom_contour2(..., gap = gap),
           geom_text_contour(..., rotate = rotate))
   } else {
      geom_contour2(...)
   }
}

# For vertical cross-sections
coord_latlev <- function(ratio = 20, ...) coord_fixed(ratio = ratio, ...)
coord_lonlev <- function(ratio = 20*4, ...) coord_fixed(ratio = ratio, ...)

lev.breaks <- c(10, 30, 100, 200, 500, 1000)
```

```{r read-ncep}
ncep.f <- memoise::memoise(function(lat = -90:40, 
                                    lon = 0:360,
                                    time = lubridate::as_datetime(c("1979-12-01", "2015-12-01")),
                                    level = 10:1000,
                                    vars = "gh"){
   subset <- list(lat = lat, lon = lon, level = level,
                  time = time)
   n <- ReadNetCDF("DATA/hgt.mon.mean.nc", vars = c(gh = "hgt"),
                   subset = subset) %>% 
      setnames(., c("level", "time"), c("lev", "date"))
   if ("u" %in% vars) {
      n[, u := ReadNetCDF("DATA/uwnd.mon.mean.nc", out = "vector",
                          subset = subset)[[1]]] 
   }
   if ("v" %in% vars) {
      n[, v := ReadNetCDF("DATA/vwnd.mon.mean.nc", out = "vector",
                          subset = subset)[[1]]]
   }
   n[, date := as.Date(date[1]), by = date]
   return(n)
}, cache = memoise::cache_filesystem(".rcache"))
ncep <- ncep.f()
```

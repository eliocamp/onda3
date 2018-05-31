library(data.table)
library(magrittr)

# http://nsidc.org/data/NSIDC-0116
setwd("DATA/ice_drift")
map <- fread("south_x_y_lat_lon.txt")
colnames(map) <- c("x", "y", "lat", "lon")

files <- list.files("bin")
years <- stringr::str_sub(files, 22, 25)
months <- stringr::str_sub(files, 27, 28)

dates <- as.Date(paste0(years, "-", months, "-01"))
grid <- setDT(expand.grid(x = 0:320, y = 0:320, 
                          time = dates[1],
                          u = 0,
                          v = 0,
                          q = 0))
grid <- grid[map, on = c("x", "y")]
setcolorder(grid, c(1:2, 7:8, 3:6))

ice.drift <- lapply(seq_along(files), function(i) {
   g <- readBin(file(paste0("bin/", files[i]), "rb"), "int", n = 800000, 
                size = 2, endian = "little") 
   grid <- grid[, `:=`(time = dates[i],
                       u = metR::JumpBy(g, 3), 
                       v = metR::JumpBy(g, 3, start = 2),
                       q = metR::JumpBy(g, 3, start = 3))]
 
   closeAllConnections()
   copy(grid[q > 0])
})
ice.drift <- rbindlist(ice.drift)

saveRDS(ice.drift, file = "../ice.drift.rda")


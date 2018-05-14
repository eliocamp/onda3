library(data.table)

ReadSeaIceBIN <- function(file, out = c("data", "land")) {
   con = file(file, "rb")
   header <- as.list(readBin(con, "int", size = 1, n = 300))
   
   # names(header) <- c("na.value", "ncol", "nrow", "", "lat.max", "orientation",
   #                    "", "j0", "i0", "instrument", "data", "julian.day.start", "hour.start" ,
   #                    "minute.start", "julian.day.end", "hour.end", "minute.end", "year",
   #                    "julian.day", "channel", "scaling.factor")
   # 
   # # other <- as.list(readBin(file, "single", size = 6, n = 3))
   # filename <- readBin(file, "character", size = 21, n = 1)
   # title <- readBin(file, "character", size = 80, n = 1)
   # other <- readBin(file, "character", size = 70, n = 1)
   
   data <- as.numeric(readBin(con, "int", size = 1, signed = FALSE,
                              n = 316*332))
   close(con)
   # data[data < 0] <- data[data < 0] + 128
   dim(data) <- c(316, 332)
   dimnames(data) <- list(x = 1:316, y = 1:332)
   
   data <- as.data.table(data.table::melt(data, value.name = "concentration"))
   
   if (out[1] == "data") {
      # ice <- data[value %between% c(0, 250) | value == 255][value == as.numeric(header$na.value), value := NA]
      ice <- copy(data)[!(concentration %between% c(0, 250)), concentration := NA]
      # Scale data
      ice[, concentration := concentration/250]
      date <- stringr::str_sub(basename(file), 4, 9)
      year <- stringr::str_sub(date, 1, 4)
      month <- stringr::str_sub(date, 5, 6)
      
      ice[, date := lubridate::ymd(paste0(year, "-", month, "-", 01))]
      
   } else {
      land.mask <- suppressWarnings(data[concentration == 254][, concentration := TRUE])
   }  
}


sea <- rbindlist(lapply(list.files("DATA/seaice", full.names = TRUE), ReadSeaIceBIN))

file <- "DATA/seaice/nt_198601_n07_v1.1_s.bin"
land.mask <- ReadSeaIceBIN(file, out = "land")

sea[, y := -y + max(y)]

sea[, x1 := - 3950*1000 + 2*3950*1000*x/316, by = x]
sea[, y1 := - 3950*1000 + (4350+3950)*1000*y/332, by = y]
proj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"
coord <- sea[date == date[1]] %>% 
   .[, c("lon", "lat") := proj4::project(list(x1, y1), proj = proj, 
                                         inverse = TRUE)] %>% 
   .[, .(x, y, lon, lat)]

sea[, c("lon", "lat") := coord[, .(ConvertLongitude(lon), lat)]]


saveRDS(sea, "DATA/seaice.Rds")
saveRDS(land.mask, "DATA/landmask.Rds")

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
      
      ice[, time := lubridate::ymd(paste0(year, "-", month, "-", 01))]
      
   } else {
      land.mask <- suppressWarnings(data[concentration == 254][, concentration := TRUE])
   }  
}

xy2lonlat <- function(x, y) {
   data <- data.table(x, y)
   
   data[, y1 := as.numeric(-y + max(y))]
   data[, y1 := - 3950*1000 + 25000*y1]
   data[, x1 := - 3950*1000 + 25000*x]
   
   datau <- unique(data)
   
   proj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"
   datau[, c("lon", "lat") := proj4::project(list(x1, y1), proj = proj, 
                                             inverse = TRUE)]
   datau[, lon := metR::ConvertLongitude(lon, from = 180)]
   
   as.list(datau[data, on = c("x", "y")][, .(lon, lat)])
}


sea <- rbindlist(lapply(list.files("DATA/seaice", full.names = TRUE), ReadSeaIceBIN))

sea[, conc := concentration]
sea[is.na(conc), conc := 0]
sea[, keep := !all(conc == 0), by = .(x, y)]
sea <- sea[keep == TRUE, .(time, x, y, concentration)]

file <- "DATA/seaice/nt_198601_n07_v1.1_s.bin"
land.mask <- ReadSeaIceBIN(file, out = "land")

sea[, c("lon", "lat") := xy2lonlat(x, y)]

saveRDS(sea, "DATA/seaice.Rds")
saveRDS(land.mask, "DATA/landmask.Rds")

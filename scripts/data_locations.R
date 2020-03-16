ERAI <- function() {
   file <- here::here("DATA", "reanalysis", "ERA-Interim", "erai.mon.mean.nc")
   checkmate::assert_access(file, access = "r")
   file
}

ERA20 <- function() {
   file <- here::here("DATA", "ERA-20C", "era20c.mon.mean.nc")
   checkmate::assert_access(file, access = "r")
   file
}

ERA5 <- function(temporal = c("mon", "day"), vertical = c("pl", "sl"), var) {
   vertical <- vertical[1]
   checkmate::assert_choice(vertical, c("pl", "sl"))
   temporal <- temporal[1]
   checkmate::assert_choice(temporal, c("mon", "day"))
   if (temporal == "mon") {
      if (vertical[1] == "pl") {
         file <- here::here("DATA", "reanalysis", "ERA5", "mon", "era5.mon.mean.nc")
      } else {
         file <- here::here("DATA", "reanalysis", "ERA5", "mon", "era5sl.mon.mean.nc")
      }   
   } else {
      stop("not implemented")
   }
   
   checkmate::assert_access(file, access = "r")
   file
}


NCEP <- function(temporal = "mon", vertical = c("pl", "flux", "sfc", "sigma")) {
   vertical <- vertical[1]
   checkmate::assert_choice(vertical, c("pl", "flux", "sfc", "sigma"))
   temporal <- temporal[1]
   checkmate::assert_choice(temporal, "mon")
   
   if (vertical == "pl") {
      file <- here::here("DATA", "reanalysis", "NCEP", temporal, "ncep.mon.mean.nc")
   } else {
      file <- here::here("DATA", "reanalysis", "NCEP", temporal, paste0("ncep-", vertical, ".mon.mean.nc"))
   }
   
   checkmate::assert_access(file, access = "r")
   file
}

OLR <- function() {
   file <- here::here("DATA", "reanalysis", "NOAA", "olr.mon.mean.nc")
   checkmate::assert_access(file, access = "r")
   file
}

CMAP <- function() {
   file <- here::here("DATA", "reanalysis", "CMAP", "precip.mon.mean.nc")
   checkmate::assert_access(file, access = "r")
   file
}

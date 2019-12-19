


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

ERA5 <- function(which = c("pl", "sl")) {
   if (which[1] == "pl") {
      file <- here::here("DATA", "reanalysis", "ERA5", "era5.mon.mean.nc")
   } else {
      file <- here::here("DATA", "reanalysis", "ERA5", "era5sl.mon.mean.nc")
   }
   checkmate::assert_access(file, access = "r")
   file
}


NCEP <- function(which = c("pl", "flux", "sfc")) {
   which <- which[1]
   if (which == "pl") {
      file <- here::here("DATA", "reanalysis", "NCEP", "ncep.mon.mean.nc")
   } else {
      file <- here::here("DATA", "reanalysis", "NCEP", paste0("ncep-", which, ".mon.mean.nc"))
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

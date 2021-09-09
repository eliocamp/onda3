print.nc_file <- function(x, ...) {
   print(metR::GlanceNetCDF(x))
}

SST <- function() {
   file <- here::here("DATA", "reanalysis", "ERSST_V5", "sst.ersst.mnmean.nc")
   checkmate::assert_access(file, access = "r")
   class(file) <- c("nc_file", class(file))
   return(file)
}

ERAI <- function() {
   file <- here::here("DATA", "reanalysis", "ERA-Interim", "erai.mon.mean.nc")
   checkmate::assert_access(file, access = "r")
   class(file) <- c("nc_file", class(file))
   return(file)
}

ERA20 <- function() {
   file <- here::here("DATA", "reanalysis", "ERA20C", "era20c.mon.mean.nc")
   checkmate::assert_access(file, access = "r")
   class(file) <- c("nc_file", class(file))
   return(file)
}

available_era5 <- list.files(here::here("DATA", "reanalysis", "ERA5", "day"))
vars <- gsub(".day.mean.nc", "", gsub("era5.", "", available_era5))

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
   } else if (temporal == "day") {
      available_era5 <- list.files(here::here("DATA", "reanalysis", "ERA5", "day"))
      vars <- gsub(".day.mean.nc", "", gsub("era5.", "", available_era5))
      var <- var[1]
      checkmate::assert_choice(var, vars)
      
      file <- paste0("era5.", var[1], ".day.mean.nc")
      file <- here::here("DATA", "reanalysis", "ERA5", "day", file)
   }
   
   checkmate::assert_access(file, access = "r")
   
   class(file) <- c("nc_file", class(file))
   return(file)
}

formals(ERA5)$var <- vars
remove(vars)

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
   
   class(file) <- c("nc_file", class(file))
   return(file)
}

OLR <- function() {
   file <- here::here("DATA", "reanalysis", "NOAA", "olr.mon.mean.nc")
   checkmate::assert_access(file, access = "r")
   class(file) <- c("nc_file", class(file))
   return(file)
}

CMAP <- function() {
   file <- here::here("DATA", "reanalysis", "CMAP", "precip.mon.mean.nc")
   checkmate::assert_access(file, access = "r")
   class(file) <- c("nc_file", class(file))
   return(file)
}

PSA <- function() {
   file <- here::here("DATA", "psa-pattern.csv")
   checkmate::assert_access(file, access = "r")
   return(data.table::fread(file))
}


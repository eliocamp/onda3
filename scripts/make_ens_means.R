library(data.table)
library(magrittr)
# "Sincroniza" los datos de pikachu -------------------------------------------

base_dir <- "/home/elio/CMIP6"
dir_to <- "/home/elio/DATOS/CMIP6"

mount <- mountr::sshfs_mount("elio.campitelli", 
                             "pikachu.cima.fcen.uba.ar", 
                             "/datos3/CMIP6",
                             base_dir)

cimadata::cmip_available(base_dir) %>% 
   as.data.table() %>% 
   .[experiment_id == "historical" & variable_id == "zg"] %>% 
   .[, file_to := gsub(base_dir, dir_to, file)] %>% 
   .[, if (!file.exists(file_to)) file.copy(file, file_to), by = file]


mountr::sshfs_unmount(mount)

# Hace la media del ensamble y medias estacionales  ----------------------------
files_in <- cimadata::cmip_available(dir_to) %>% 
   as.data.table() %>% 
   .[experiment_id == "historical" & variable_id == "zg"] %>% 
   .$file

files_out <- paste0("DATA/", basename(files_in))


make_mean <- function(file_in, file_out, ensemble_mean = TRUE, verbose = FALSE) {
   library(reticulate)
   py_install("netcdf4")
   xr <- import("xarray")
   
   if (file.exists(file_out)) {
      if (verbose) cimadata:::message_time("Skipping ", file_in)
      return(file_out)
   }
   
   
   if (verbose) cimadata:::message_time("Processing ", file_in)
   
   if (verbose) cimadata:::message_time("Reading data")
   dt <- xr$open_dataset(file_in)$
      sel(plev = c(50*10, 200*10))$
      zg
   
   if (verbose) cimadata:::message_time("Downsampling to seasonal")
   dt <- dt$resample(time = "QS")$
      mean(dim = "time")
   
   if (ensamble_mean) {
      if (verbose) cimadata:::message_time("Averaging ensemble")      
      dt <- dt$mean(dim = "ensemble")
   }
   
   if (verbose) cimadata:::message_time("Saving to ", file_out)      
   dt$to_netcdf(file_out)
   
   return(file_out)
}

files <- purrr::map2(files_in, files_out, make_mean, verbose = TRUE)


xarray_to_dt <- function(xarray) {
   dims <- lapply(xarray$coords$dims, function(d) xarray$coords$get(d)$values)
   names(dims) <- xarray$coords$dims
   dt <- xarray$values 
   dimnames(dt) <- dims
   
   dt <- data.table::setDT(reshape2::melt(dt))
   
   if ("time" %in% colnames(dt))  {
      dt[, time := lubridate::as_datetime(time)][]  
   }
   
   dt
}


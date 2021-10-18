library(rnoaa)
library(data.table)
library(magrittr)

stations <- ghcnd_stations() %>% 
   as.data.table() %>% 
   .[latitude <= -20 & element %in% c("TAVG", "PRCP") & last_year >= 2018 & first_year <= 1980] %>% 
   .[, .(id, lat = latitude, lon = longitude, variable = element)] %>% 
   unique()

station_data <- ghcnd_search(unique(stations$id), var = c("tavg", "prcp")) %>% 
   rbindlist(idcol = "variable") %>% 
   as.data.table() 

station_data_clean <- station_data[unique(stations[, .(id, lon, lat)]), on = "id"] %>% 
   .[, time := lubridate::as_datetime(date[1]), by = date] %>% 
   .[, date := NULL] %>% 
   .[, lon := metR::ConvertLongitude(lon)] %>% 
   setnames("tavg", "value") %>% 
   .[variable != ""] %>% 
   .[, value := as.double(value)] %>% 
   .[variable == "prcp", value := value/10] %>% 
   .[variable == "tavg", value := value/10] %>% 
   .[qflag == " ", value_clean := value] 

fwrite(station_data_clean, here::here("DATA", "ghcnd_stations.csv"), yaml = TRUE)

      
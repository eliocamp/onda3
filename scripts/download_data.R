setwd("~/Documents/CONICET/onda3/")

# Ice Concentration
f <- readLines("DATA/ice_files.txt")
urls <- paste0('https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0051_gsfc_nasateam_seaice/final-gsfc/south/monthly/',
                           f)

r <- lapply(seq_along(urls), function(x) {
   command <- paste0('wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off ',
                     urls[x])
   system(command)
})

setwd("DATA/ice_drift")
# Ice drift
years <- 1980:2017
months <- 1:12
comb <- as.data.table(expand.grid(month = months, year = years))
comb[, month := ifelse(month < 10, paste0("0", month), month)]
base <- "https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0116_icemotion_vectors_v3/data/south/means/"
r <- lapply(seq_len(nrow(comb)), function(i) {
   year <- comb$year[i]
   month <- comb$month[i]
   file <- paste0("icemotion.grid.month.", year, ".", month, ".s.v3.bin")
   url <- paste0(base, year, "/", file)
   command <- paste0('wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off ',
                     url)
   system(command)
})

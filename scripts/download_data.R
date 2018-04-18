setwd("~/Documents/CONICET/onda3/")

f <- readLines("DATA/ice_files.txt")

urls <- paste0('https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0051_gsfc_nasateam_seaice/final-gsfc/south/monthly/',
                           f)

r <- lapply(seq_along(urls), function(x) {
   command <- paste0('wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off ',
                     urls[x])
   system(command)
})



library(magrittr)

file <- "23-mas_SON.Rmd"
seasons <- c("DJF", "MAM", "JJA", "SON")
PCs <- 1:2

expand.grid(season = seasons, PC = PCs) %>% 
   .[1:2, ] %>% 
   purrr::pmap(~ rmarkdown::render(file, 
                                   params = list(season = .x, PC = .y), 
                                   output_file = paste0("24 - ", .x, " ", .y)
                                   ))

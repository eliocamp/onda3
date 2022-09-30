

trimestres <- c("SON", "MAM", "JJA", "DJF")


for (trimestre in trimestres) {
   message("Renderizando ", trimestre)
   rmarkdown::render("35-cEOF-CMIP6-superensemble.Rmd", 
                     output_file = paste0("35-cEOF-CMIP6-superensemble-", trimestre),
                     params = list(season = trimestre))
}

message("Listo!")


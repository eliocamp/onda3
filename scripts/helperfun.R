## Archivo con funciones

#Librerías
library(ggplot2)
library(ggforce)
library(stringi)
library(ggthemes)
library(magrittr)
library(data.table)
library(viridis)
library(lubridate)
library(directlabels)
library(akima)
library(compiler)
library(RColorBrewer)
enableJIT(0)

# Mapa
BuildMap <- function(res = 1, smooth = 1, pm = 180,
                     countries = FALSE, ...) {
   # Caraga datos de mapas y cambia el meridiano principal, suaviza y cambia
   # resolución para plotear más rápido.
   # Entra:
   #   res: resolución (cuantos puntos se eliminan entre los que quedan)
   #   smooth: factor de suavizado (ventana del promedio corrido)
   #   pm: longitud del meridiano central
   #   countries: ¿datos a nivel país?
   # Sale:
   #   un data.table con las coordenadas de cada polígono y su grupo
   library(data.table)
   if (countries) {
      library(rworldxtra)
      data("countriesLow")
      m <- countriesLow
   } else {
      library(rworldmap)
      data(coastsCoarse)
      m <- coastsCoarse
   }
   
   m <- as.data.table(fortify(m))
   m[, group := as.numeric(group)]
   m[, id := as.numeric(id)]
   
   # Cambio el prime meridian.
   m2 <- copy(m)
   m2[, long := long + 360]
   m2[, group := group + max(group) + 1]
   m <- rbind(m, m2)
   m <- m[long >= pm - 180 & long <= 180 + pm]
   
   m[, MO := max(order), by = group]
   
   # Suavizo.
   if (smooth > 1) {
      cut <- max(smooth, res)
      notsmooth <- m[, .N, by = group][N < cut + 4, group]
      m <- m[!(group %in% notsmooth)]    # saco los grupos muy chicos
      
      m[order != 1 & order != MO,
        `:=`(long = zoo::rollmean(long, smooth, fill = "extend"),
             lat = zoo::rollmean(lat, smooth, fill = "extend")),
        by = group]
   }
   # Bajo la resolución.
   suppressWarnings(m[, keep := c(1, rep(0, res - 1)), by = group])
   m <- m[order == 1 | order == MO | keep == 1, .SD, by = group]
   
   return(m)
}

geom_map2 <- function(map, size = 0.2, color = "gray50") {
   # Un geom_path con defaults copados para agregar mapas
   g <- geom_path(data = map, aes(long, lat, group = group),
                  inherit.aes = F, color = color, size = size)
   return(g)
}


# Repetir vector en negativo
fill_neg <- function(vector) {
   c(-vector[order(-vector)], vector)
}



# Nombres de los meses en español
month.abb_sp <- c("Ene", "Feb", "Mar", "Abr", "May", "Jun",
                  "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")
names(month.abb_sp) <- as.character(1:12)

names(month.abb) <- as.character(1:12)


# Para interpolación en data table
Interpolate.DT <- function(z, x, y, yo = unique(y), xo = unique(x), ...){
   na <- is.na(z)
   int <- akima::interp(x = x[!na], y = y[!na], z = z[!na], yo = yo, xo = xo, ...)
   names <- c(deparse(substitute(x)),
              deparse(substitute(y)),
              deparse(substitute(z)))    # muy feo, sí
   r <- with(int, {
      grid <- expand.grid(x, y)
      r <- list(grid[,1], grid[, 2], c(z))
      names(r) <- names
      return(r)
   })
}

# Función que hace autocorrelograma y su test según Anderson o large lag.
acf.sig <- function(x, lag.max=0.3*length(x), alpha = 0.05,
                    method=c("anderson","large.lag", "salas"), sided="one") {
   autocor <- acf(x, lag.max=lag.max, plot = F)$acf
   N <- length(x)
   e <- -1/(N-1)
   if (method[1]=="anderson"){
      var <- (N-2)/(N-1)^2
   } else if (method[1]=="large.lag"){
      var <- vector()
      for (i in 1:length(autocor)){
         v <- ifelse(i==1,1/N, 1/N * (1+2*sum(autocor[1:i-1]^2)))
         var<- c(var, v)
      }
   } else if (method[1]=="salas"){
      var <- vector()
      
      for (i in 1:length(autocor)){
         v <- (N-1-i)/(N-i)^2
         var<- c(var, v)
         e <- -1/(N-i)
      }
   }
   if (sided=="one"){
      a <- alpha
      q <- qnorm(a, lower.tail=F)
      sigupp <- e + sqrt(var)*q
      ret <- data.table(lag=0:lag.max, acf=autocor, sig.cut=sigupp)
   } else if (sided == "two"){
      a <- alpha/2
      q <- qnorm(a, lower.tail=F)
      sigupp <- e+sqrt(var)*q
      siginf <- e-sqrt(var)*q
      ret <- data.table(lag=0:lag.max, acf=autocor, upp.sig.cut=sigupp, low.sig.cut=siginf)
   }
   ret
}



# Convierte la salida de la función fft en un formato
# legible por humanos.

convert.fft <- function(cs, sample.rate=1, full=T) {
   distance.center <- function(c) Mod(c)
   angle <- function(c) Arg(c)
   is.even <- function(x) ceiling(x/2) == x/2
   N <- length(cs)
   if (full==T){
      nyq <- ifelse(is.even(N), N/2+1, (N+1)/2)
      cs <- cs[2:nyq]
   }
   NP <- length(cs)
   cs <- cs / N # normalize
   
   df <- data.frame(cycle    = 1:(NP),
                    freq     = 1:(NP) * sample.rate / N,
                    per      = N/(1:(NP) * sample.rate),
                    ampl = sapply(cs, Mod),
                    delay    = sapply(cs, angle),
                    spect    = sqrt(sapply(cs, Mod)),
                    comp = cs)
   
   non.unique <- ifelse(is.even(NP), NP-1, NP)
   df$ampl [1:non.unique] <- df$ampl[1:non.unique]*2
   df
}

# ¿por qué está esto acá? Es una función que convierte una columna con n factores
# en n columnas con 0 o 1. ¿No es un dcast? Que la pedo que estoy...
factor2cols <- function(x, column, factors) {
   column <- deparse(substitute(column))
   for(h in factors){
      x[, (h) := ifelse(get(column) == h, 1, 0)]
   }
   x[is.na(get(column)), (factors) := 0]
}

# convierte una fecha en formato AAAA-MM-DD en mes con factor y ordenado según estaciones
date_month2factor <- function(x) {
   factor(as.numeric(stringi::stri_sub(x, 6, 7)), levels = c(12, 1:11), ordered = T)
}

#### Esto es re útil y lo tengo que poner en metR. 
#### Técnicametne es igual que broom::tidy, pero ese parece que no funciona bien 
#### con RcppArmadillo::fastLm
ExtractLm <- function(model) {
   # Extrae el estimador y el error estándar de cada regressor para un modelo
   # lineal (¿o cualquier modelo?).
   # Entra:
   #   model: un modelo lineal (o no)
   # Sale:
   #   una lista con 3 vectores: el nombre de los elementos, los estimadores y
   #   el error estándar.
   a <- summary(model)
   return(list(regressor = rownames(a$coefficients),
               estimate  = unname(a$coefficients[, 1]),
               se        = unname(a$coefficients[, 2])))
}

WaveFlux <- function(psi, p = 250, a = 6371000) {
   k <- p*100/(a^2*2000)
   psi <- copy(psi)
   psi[, c("psi.dlon", "psi.dlat") := Derivate(psi.z ~ lon + lat,
                                               cyclical = c(TRUE, FALSE))] %>%
      .[, psi.ddlon := Derivate(psi.z ~ lon, cyclical = TRUE, order = 2),
        by = lat] %>%
      .[, psi.dlondlat := Derivate(psi.dlon ~ lat),
        by = lon] %>%
      .[, `:=`(f.lon = k/cos(lat*pi/180)*(psi.dlon^2 - psi.z*psi.ddlon),
               f.lat = k*(psi.dlon*psi.dlat - psi.z*psi.dlondlat))]
   list(f.lon = psi$f.lon, f.lat = psi$f.lat)
}

guide_colorstrip_bottom <- function(width = 25, height = 0.5, ...) {
   guide_colorstrip(title.position = "top", title.hjust = 0.5,
                    barheight = height,
                    barwidth = width, ...)
}

scale_x_longitude <- function(ticks = 60, ...) {
   metR::scale_x_longitude(ticks = ticks, breaks = seq(-180, 360, by = ticks), ...)
}
scale_s_map <- function(ylim = c(-90, 40),
                        xlim = c(0, 360)) list(scale_y_latitude(limits = ylim),
                                                     scale_x_longitude(limits = xlim))

ggplot <- function(...) {
   ggplot2::ggplot(...) + scale_linetype(guide = "none")
}

AddSuffix <- function(suffix = "") {
   function(string) {
      paste0(string, suffix)
   }
}

AddPreffix <- function(preffix = "") {
   function(string) {
      paste0(preffix, string)
   }
}

lev.lab <- AddSuffix(" hPa")
qs.lab <- AddPreffix("QS ")


yearmonth <- function(date, day = 1) {
   months <- lubridate::month(date)
   years <- lubridate::year(date)
   lubridate::ymd(paste(years, months, day, sep = "-"))
}

yearly <- function(date, day = 182) {
   years <- lubridate::year(date)
   d <- lubridate::ymd(paste(years, "01", "01", sep = "-"))
   yday(d) <- day
   d
}

geom_index.region <- function(data) {
   geom_rect(data = data, aes(xmin = latmin, xmax = latmax,
                              ymin = levmin, ymax = levmax),
             inherit.aes = F, linetype = 3, color = "black", fill = NA)
}


PeriodicWavelet <- function(x, k) {
   period <- length(x)/k
   x1 <- rep(x, 3)
   keep <- (length(x)+1):(2*length(x))
   res <- list()
   for (p in seq_along(period)) {
      w <- WaveletComp::WaveletTransform(x1, dt = 1, upperPeriod = period[p],
                                         lowerPeriod = period[p])
      res[[paste0("k", ".", k[p])]] <- w$Ampl[keep]*sd(x)
      
   }
   return(res)
}

ReconstructWavelet <- function(x, k) {
   period <- length(x)/k
   x1 <- rep(x, 3)
   keep <- (length(x)+1):(2*length(x))
   w <- WaveletComp::analyze.wavelet(data.frame(x1), make.pval = F,
                                     loess.span = 0, verbose = F)
   r <- WaveletComp::reconstruct(w, sel.period = period,
                                 plot.rec = F, verbose = F)$series$x1.r
   r[keep]
}


greater <- function(x, N) {
   r <- frank(-x, ties.method = "first")
   r <= N
}


decade <- function(year) {
   substr(year, 3, 4)
}

fsign <- function(x) {
   f <- sign(x)
   factor(ifelse(f == 0, NA, f))
}



ifelse2 <- function(x, expression, yes = NA, no = x) {
   e <- eval(parse(text = deparse(substitute(expression))), envir = environment())
   yes <- eval(parse(text = deparse(substitute(yes))), envir = environment())
   no <- eval(parse(text = deparse(substitute(no))), envir = environment())
   ifelse(e, yes, no)
}

geom_label_contour2 <- function(...) {
   list(geom_label_contour(fill = "white", label.r = unit(0, "lines"),
                           label.padding = unit(0.06, "lines"), color = NA, ...),
        geom_text_contour(..., rotate = FALSE))
}


cache.file <- function(file, expression, verbose = interactive()) {
   
   if (file.exists(file)) {
      if (verbose) message("Reading data from file.")
      return(readRDS(file))
   } else {
      if (verbose) message("Evaluating expression.")
      r <- eval(expression)
      if (verbose) message("Saving data to file.")
      saveRDS(r, file = file)
      return(r)
   }
}

mode.circular <- function(x, limits = c(0, 2/3*pi)) {
   if (length(x) > 1) {
      x1 <- c(x - limits[2], x, x + limits[2])
      keep <- (length(x) + 1):(2*length(x))
      d <- density(x1)
      y <- d$y[d$x %b% limits]
      x2 <- d$x[d$x %b% limits]
      x2[which.max(y)]
   } else {
      x
   }
}

# "Estaciones" en base a la amplitud y fase de la onda 3.
qs.season <- function(month) {
   if (metR:::.is.somedate(month)) month <- lubridate::month(month)
   
   qs.seasons <- factor(c(rep("DJFM", 3),
                          rep("A", 1),
                          rep("MJJ", 3),
                          rep("ASO", 3),
                          rep("N", 1),
                          rep("DJFM", 1)))
   
   return(factor(qs.seasons[month], levels = c("DJFM", "A", "MJJ",
                                               "ASO", "N")))
}


# Trimestres
qs.trim <- function(month) {
   if (metR:::.is.somedate(month)) month <- lubridate::month(month)
   
   qs.seasons <- factor(c(rep("JFM", 3),
                          rep(NA, 1),
                          rep("MJJ", 3),
                          rep("ASO", 3),
                          rep(NA, 2)))
   
   return(factor(qs.seasons[month], levels = c("JFM", "MJJ",
                                               "ASO")))
}



geom_contour_back <- function(..., color = "black", size = 0.2, alpha = 0.5) {
   geom_contour(..., color = color, size = size, alpha = alpha)
}
geom_label_contour_back <- function(...) {
   geom_label_contour2(..., alpha = 0.5, size = 3)
}

geom_contour_fine <- function(...) geom_contour(..., size = 0.4)
stat_contour_fine <- function(...) stat_contour(..., size = 0.4)

geom_cross <- function(x = 0, y = 0, ...) {
   list(geom_vline(xintercept = x, ...),
        geom_hline(yintercept = y, ...))
}

labeller.date <- function(sep = " - ") {
   function(s) {
      s <- as.Date(s)
      m <- month(s)
      y <- year(s)
      paste0(month.abb_sp[m], sep, y)
   }
}

no.zero_ <- function(x) {
   if (x == 0) return(".0")
   if (abs(x) < 1) {
      s <- ifelse(x < 0, "-", "")
      paste0(s, substr(abs(x), 2, nchar(x)))
   } else {
      x
   }
}

no.zero <- function(x) {
   sapply(seq_along(x), function(i) no.zero_(x[i]))
}


StatContour3 <- ggplot2::ggproto("StatContour3", metR:::StatContour2,
                                 default_aes = aes(order = calc(level), linetype = factor(-sign(calc(level))))
)

geom_contour3 <- function(mapping = NULL, data = NULL,
                          stat = "contour3", position = "identity",
                          ...,
                          color = "black",
                          size = 0.3,
                          lineend = "butt",
                          linejoin = "round",
                          linemitre = 1,
                          breaks = scales::fullseq,
                          bins = NULL,
                          binwidth = NULL,
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE) {
   layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomContour2,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
         lineend = lineend,
         linejoin = linejoin,
         linemitre = linemitre,
         breaks = breaks,
         bins = bins,
         binwidth = binwidth,
         na.rm = na.rm,
         color = color,
         size = size,
         ...
      )
   )
}


FilterWave <- function(x, k) {
   f <- fft(x)
   # Need to remove the k+1 spots (because index 1 is k = 0) 
   # and the N - k + 1 because of symmetry.
   k <- c(k + 1, length(x) - k[k != 0] + 1) 
   f[k] <- 0 + 0i
   Re(fft(f, inverse = T))/length(x)
}


FitLm <- function(y, ..., se = FALSE) {
   X <- cbind(mean = 1, ...)
   regressor <- dimnames(X)[[2]]
   a <- .lm.fit(X, y)
   estimate <- a$coefficients
   if (se == TRUE) {
   sigma <- sum(a$residuals^2)/(nrow(X) - ncol(X))
   se <- sqrt(diag(solve(t(X)%*%X)*sigma))
   return(list(regressor = dimnames(X)[[2]],
               estimate = estimate,
               se = se))
   } else {
      return(list(regressor = dimnames(X)[[2]],
                  estimate = estimate))
   }
}

CutEOF <- function(eof, pc) {
   if (is.numeric(pc)) pc <- paste0("PC", pc)
   lapply(eof, function(x) {
      x[PC %in% pc]
   })
}

PermTest <- function(y, ..., N = 10) {
   original <- FitLm(y, ..., se = FALSE)
   regressor <- original$regressor
   estimate <- original$estimate
   f <- rep(0, length(estimate))
   n <- length(y)
   p <- seq_len(n)
   set.seed(42)
   for (i in seq_len(N)) {
      y <- y[sample(p, n, replace = FALSE)]
      e <- FitLm(y, ..., se = FALSE)$estimate
      f <- f + as.numeric(abs(e) >= abs(estimate))
   }
   f <- f/N
   return(append(original, list(p.value = f)))
}


Jump <- function(x, by = 1) {
   keep <- JumpBy(unique(x), by = by)
   x[!(x %in% keep)] <- NA
   x
}

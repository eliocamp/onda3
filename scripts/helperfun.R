## Archivo con funciones

#Librerías
library(ggplot2)
library(ggforce)
library(stringi)
# library(ggthemes)
library(magrittr)
library(data.table)
library(lubridate)
library(akima)
library(compiler)
# library(RColorBrewer)
enableJIT(0)

# Mapa
BuildMap <- function(res = 1, smooth = 1, pm = 180,
                     countries = FALSE, ...) {
   # Caraga datos de mapas y cambia el meridiano principal, suaviza y cambia
   # resolución para plotear más rápido.
   # Entra:
   #   res: resolución (cuantos puntos se eliminan entre los que quedan)
   #   smooth: factor de suavifzado (ventana del promedio corrido)
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


# # Para interpolación en data table
# Interpolate.DT <- function(z, x, y, yo = unique(y), xo = unique(x), ...){
#    na <- is.na(z)
#    int <- akima::interp(x = x[!na], y = y[!na], z = z[!na], yo = yo, xo = xo, ...)
#    names <- c(deparse(substitute(x)),
#               deparse(substitute(y)),
#               deparse(substitute(z)))    # muy feo, sí
#    r <- with(int, {
#       grid <- expand.grid(x, y)
#       r <- list(grid[,1], grid[, 2], c(z))
#       names(r) <- names
#       return(r)
#    })
# }

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

scale_x_longitude <- function(ticks = 60, xwrap = c(0, 360), ...) {
   b <- seq(min(xwrap), max(xwrap), by = ticks)
   metR::scale_x_longitude(ticks = ticks, breaks = b, labels = LonLabel(b),...)
}
scale_s_map <- function(ylim = c(-90, -15), xlim = c(0, 360)) {
   list(scale_y_latitude(limits = ylim),
        scale_x_longitude(xwrap = xlim)) 
} 

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
qs.sem <- function(month) {
   if (metR:::.is.somedate(month)) month <- lubridate::month(month)
   
   qs.seasons <- factor(c(rep("DJFMA", 4),
                          rep("MJJASON", 7),
                          rep("DJFMA", 1)))
   
   return(factor(qs.seasons[month], levels = c("DJFMA", "MJJASON")))
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
                          rep("ND", 2)))
   
   return(factor(qs.seasons[month], levels = c("JFM", "MJJ",
                                               "ASO", "ND")))
}
qs.sem <- qs.trim


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
                                 default_aes = aes(order = stat(level), linetype = factor(-sign(stat(level))))
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


## Methos for EOF

cut.eof <- function(eof, n) {
   var <- attr(eof, "suffix")
   value.var <- attr(eof, "value.var")
   return(structure(lapply(as.list(eof), function(x) {
      x[as.numeric(get(var)) %in% n]
   }),
   class = c("eof", "list"),
   suffix = var,
   value.var = value.var))
}

screeplot.eof <- function(eof, n = "all") {
   var <- attr(eof, "suffix")
   r2 <- "r2"
   if (n[1] == "all") n <- as.numeric(unique(eof$sdev[[var]]))
   ggplot(eof$sdev[as.numeric(get(var)) %in% seq_along(n)], aes_(as.name(var), as.name(r2))) +
      geom_point()
}

autoplot.eof <- function(eof, n = "all") {
   screeplot(eof, n)
}

predict.eof <- function(eof, n = NULL) {
   ` %>% ` <- magrittr::`%>%`
   if (!inherits(eof, "eof")) {
      stop("eof must be an EOF object")
   }
   
   if(!is.null(n)) eof <- cut(eof, n)
   
   value.var <- attr(eof, "value.var")
   pc <- attr(eof, "suffix")
   
   right.vars <- colnames(eof$right)[!(colnames(eof$right) %in% c(pc, value.var))]
   right.formula <- as.formula(paste0(pc, " ~ ", paste0(right.vars, collapse = "+")))
   
   right <- eof$right %>% 
      .[eof$sdev, on = pc] %>% 
      .[, (value.var) := get(value.var)*sd] %>% 
      metR:::.tidy2matrix(right.formula, value.var)
   
   left.vars <- colnames(eof$left)[!(colnames(eof$left) %in% c(pc, value.var))]
   left.formula <- as.formula(paste0(pc, " ~ ", paste0(left.vars, collapse = "+")))
   left <- metR:::.tidy2matrix(eof$left, left.formula, value.var)
   
   dt <- cbind(.extend.dt(left$coldims, each = nrow(right$coldims)),
               .extend.dt(right$coldims, n = nrow(left$coldims)),
               c(t(right$matrix)%*%left$matrix))
   colnames(dt)[length(colnames(dt))] <- value.var
   return(dt)
}

.extend.dt <- function(dt, n = NULL, each = NULL) {
   if (!is.null(n)) {
      r <- as.data.table(lapply(dt, rep, n = n))
   } else {
      r <- as.data.table(lapply(dt, rep, each = each))
   }
   r
}

labeller <- function(...) {
   UseMethod("labeller")
}

labeller.default <- function(...) {
   ggplot2::labeller(...)
}

labeller.eof <- function(eof, sep = " - ") {
   var <- attr(eof, "suffix")
   with(eof$sdev, setNames(paste0(get(var), sep, scales::percent(r2)),
                           get(var)))
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

Detrend <- function(y, x) {
   nas <- is.na(y)
   m <- mean(y, na.rm = TRUE)
   if (!hasArg(x)) x <- seq_along(y)
   y[!nas] <- .lm.fit(cbind(1, x[!nas]), y[!nas])$residuals
   return(y + m)
}

Jump <- function(x, by = 1) {
   keep <- JumpBy(unique(x), by = by)
   x[!(x %in% keep)] <- NA
   x
}

shift2 <- function(x, n = 1L, fill = NA, give.names = FALSE) {
   type <- ifelse(n > 0, "lead", "lag")
   data.table::shift(x, abs(n), fill, type, give.names)
}

shiftcor <- function(x, y, lags) {
   vapply(lags, function(i) cor(x, shift2(y, i), use = "complete.obs"), 1)
}

BuildEOF <- function(formula, value.var = NULL, data = NULL, n = 1, 
                     rotate = FALSE) {
   
   if (!is.null(value.var)) {
      if (is.null(data)) stop("data must not be NULL if value.var is NULL",
                              .call = FALSE)
      data <- copy(data)
      f <- as.character(formula)
      f <- stringr::str_replace(f, "~", "\\|")
      formula <- Formula::as.Formula(paste0(value.var, " ~ ", f))
   }
   
   if (is.null(data)) {
      formula <- Formula::as.Formula(formula)
      data <- as.data.table(eval(quote(model.frame(formula, data  = data))))
   }
   
   f <- as.character(formula)
   f <- stringr::str_split(f,"~", n = 2)[[1]]
   dcast.formula <- stringr::str_squish(f[stringr::str_detect(f, "\\|")])
   dcast.formula <- as.formula(stringr::str_replace(dcast.formula, "\\|", "~"))
   
   value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])
   
   g <- metR:::.tidy2matrix(data, dcast.formula, value.var)
   
   if (is.null(n)) n <- seq_len(min(ncol(g$matrix), nrow(g$matrix)))
   
   if (requireNamespace("irlba", quietly = TRUE) &
       max(n) < 0.5 *  min(ncol(g$matrix), nrow(g$matrix))) {
      set.seed(42)
      eof <- irlba::irlba(g$matrix, nv = max(n), nu = max(n), rng = runif)
   } else {
      eof <- svd(g$matrix, nu = max(n), nv = max(n))
      eof$d <- eof$d[1:max(n)]
   }
   eof$D <- diag(eof$d, ncol = max(n), nrow = max(n))
   
   if (rotate == TRUE & max(n) > 1) {
      # Rotation
      loadings <- t(with(eof, D%*%t(v)))
      scores <- eof$u
      R <- varimax(loadings, normalize = FALSE)
      eof$u <- eof$u%*%R$rotmat
      
      # Recover rotated V and D matrixs
      loadings <- R$loadings
      class(loadings) <- "matrix"
      eof$d <- sqrt(apply(loadings, 2, function(x) sum(x^2)))
      eof$v <- t(diag(1/eof$d)%*%t(loadings))
   }
   
   c(with(eof, u%*%D%*%t(v)))
}



stat_contour4 <- function(mapping = NULL, data = NULL,
                          geom = "contour", position = "identity",
                          ...,
                          breaks = scales::fullseq,
                          bins = NULL,
                          binwidth = NULL,
                          na.rm = FALSE,
                          circular = NULL,
                          show.legend = NA,
                          inherit.aes = TRUE) {
   layer(
      data = data,
      mapping = mapping,
      stat = StatContour4,
      geom = geom,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
         na.rm = na.rm,
         breaks = breaks,
         bins = bins,
         binwidth = binwidth,
         circular = circular,
         ...
      )
   )
}

#' @rdname geom_contour2
#' @usage NULL
#' @format NULL
#' @export
StatContour4 <- ggplot2::ggproto("StatContour4", Stat,
                                 required_aes = c("x", "y", "z"),
                                 default_aes = ggplot2::aes(order = ..level..),
                                 setup_params = function(data, params) {
                                    # Check is.null(breaks) for backwards compatibility
                                    if (is.null(params$breaks)) {
                                       params$breaks <- scales::fullseq
                                    }
                                    
                                    if (is.function(params$breaks)) {
                                       # If no parameters set, use pretty bins to calculate binwidth
                                       if (is.null(params$bins) && is.null(params$binwidth)) {
                                          params$binwidth <- diff(pretty(range(data$z), 10))[1]
                                       }
                                       # If provided, use bins to calculate binwidth
                                       if (!is.null(params$bins)) {
                                          params$binwidth <- diff(range(data$z)) / params$bins
                                       }
                                       
                                       params$breaks <- params$breaks(range(data$z), params$binwidth)
                                    }
                                    return(params)
                                    
                                 },
                                 compute_group = function(data, scales, bins = NULL, binwidth = NULL,
                                                          breaks = scales::fullseq, complete = FALSE,
                                                          na.rm = FALSE, circular = NULL) {
                                    
                                    if (!is.null(circular)) {
                                       # M <- max(data[[circular]]) + resolution(data[[circular]])
                                       data <- RepeatCircular(data, circular)
                                    }
                                    
                                    contours <- as.data.table(.contour_lines(data, breaks, complete = complete))
                                    
                                    # contours <- .order_contour(contours, setDT(data))
                                    
                                    return(contours)
                                 }
)


.contour_lines <- function(data, breaks, complete = FALSE) {
   
   cl <- setDT(contoureR::getContourLines(
      x = data$x, y = data$y, z = data$z, levels = breaks))
   
   if (length(cl) == 0) {
      warning("Not possible to generate contour data", call. = FALSE)
      return(data.frame())
   }
   setnames(cl, c("z", "Group", "PID"), c("level", "group", "piece"))
   return(cl)
}


Smooth2D <- function(formula, x.out = 64, y.out = 64, data = NULL, ...) {
   dep.names <- formula.tools::lhs.vars(formula)
   if (length(dep.names) == 0) stop("LHS of formula must have at least one variable")
   
   ind.names <- formula.tools::rhs.vars(formula)
   if (length(ind.names) > 2) {
      stop("RHS of formula must be of the form x + y")
   }
   
   formula <- Formula::as.Formula(formula)
   data <- as.data.table(eval(quote(model.frame(formula, data = data,
                                                na.action = NULL))))
   
   loc <- setDF(data)[, ind.names]
   # for (v in seq_along(dep.names)) {
   value.var <- dep.names[1]
   
   sm <- fields::smooth.2d(data[[value.var]], loc, 
                           nrow = x.out, ncol = y.out, 
                           ...)
   if (!is.finite(diff(range(sm$z)))) {
      stop(paste0("smoothing failed for ", value.var, " use a bigger smooth parameter"))
   }
   
   dimnames(sm$z) <- with(sm, setNames(list(x, y), ind.names))
   
   z <- setDT(melt(sm$z, value.name = value.var))
   # set(loc, NULL, value.var, z)
   # }
   return(z)
}


seq_range <- function(x, ...) {
   r <- range(x)
   seq(r[1], r[2], ...)
}

fft2 <- function(x, k) {
   f <- fft(x)/length(x)
   f[-1] <- f[-1]*2
   return(list(R = Re(f[k + 1]), 
               I = Im(f[k + 1])))
}


StatRasa <- ggplot2::ggproto("StatRasa", Stat,
                             compute_group = function(data, scales, fun, fun.args) {
                                args <- formals(fun)
                                
                                for (i in seq_along(fun.args)) {
                                   if (names(fun.args[i]) %in% names(fun.args)) {
                                      args[[names(fun.args[i])]] <- fun.args[[i]]
                                   }
                                }
                                
                                formals(fun) <- args
                                fun(data)
                             })

stat_rasa <- function(mapping = NULL, data = NULL,
                      geom = "point", 
                      position = "identity",
                      fun = NULL,
                      ...,
                      show.legend = NA,
                      inherit.aes = TRUE) {
   if (!is.function(fun)) stop("fun must be a function")
   
   fun.args <- match.call(expand.dots = FALSE)$`...`
   layer(
      data = data,
      mapping = mapping,
      stat = StatRasa,
      geom = geom,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      check.aes = FALSE,
      check.param = FALSE,
      params = list(
         fun = fun, 
         fun.args = fun.args,
         na.rm = FALSE,
         ...
      )
   )
}


SmoothContour <- function(data, nx = 64, ny = 64, breaks, 
                          smooth = 0.05) {
   data <- data[complete.cases(data), ]
   M <- max(diff(range(data$x)), diff(range(data$y)))
   theta <- smooth*M
   sm <- setDF(Smooth2D(z ~ x + y, data = data, x.out = nx, y.out = ny,
                        theta = theta))
   if (is.function(breaks)) {
      breaks <- breaks(sm$z)      
   }
   
   contours <- as.data.table(.contour_lines(sm, breaks,
                                            complete = FALSE))
   
   if (length(contours) == 0) {
      warning("Not possible to generate contour data", call. = FALSE)
      return(data.frame())
   }
   
   contours <- metR:::.order_contour(contours, setDT(sm))
   return(contours)
}




notify <- function(title = "title", text = NULL, time = 2) {
   time <- time*1000
   system(paste0('notify-send "', title, '" "', text, '" -t ', time, ' -a rstudio'))
}

notify_after <- function(expression, ...) {
   expression <- eval(expression)
   notify(title = "Run\\ ended", ...)
   return(expression)
}

# Plus-minus functions
pm <- function(x) {
   x <- abs(x)
   c(-x, x)
}


`%pm%` <- function(x, y) {
   c(x - y, x + y)
}


theilsen <- function(formula, data, repeated = TRUE, weights = NULL) {
   mblm::mblm(formula, dataframe = data, repeated = repeated)
}

prob.t <- function(estimate, se, df) {
   pt(abs(estimate)/se, df, lower.tail = FALSE)
}

cut_obs <- function(x, n, keep = "first") {
   n_groups <- floor(length(x)/n)
   groups <- rep(seq_len(n_groups), each = n)
   nas_length <- length(x) - length(groups)
   
   if (keep == "first") {
      groups <- c(groups, rep(NA, nas_length))
   } else {
      groups <- c(rep(NA, nas_length), groups)
   }
   return(groups)
}

MakeCircle <- function(r, x0 = 0, y0 = 0, n = 40) {
   data <- data.table(r = r, x0 = x0, y0 = y0)
   theta <- seq(0, 360, length.out = n)
   data[, .(x = r*cos(theta*pi/180) + x0,
            y = r*sin(theta*pi/180) + y0), by = .(r)]
}

MakeLine <- function(angle, r = 1) {
   data <- data.table(angle = angle, x = 0, y = 0)
   data[, `:=`(xend = r[1]*cos(angle*pi/180),
               yend = r[1]*sin(angle*pi/180))]
}

angle <- function(a, b) acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )

angle2 <- function(M, N){
   atan2(N[2], N[1]) - atan2(M[2], M[1]) 
}



dervangle <- function(x, y) {
   N <- length(x)
   a_prev <- atan2(y[c(N, 1:(N-1))],  x[c(N, 1:(N-1))])
   a_next <- atan2(y[c(2:N, 1)], x[c(2:N, 1)])
   
   dxdy <- (a_next - a_prev)
   dxdy[dxdy > pi] <- dxdy[dxdy > pi] - 2*pi
   dxdy[dxdy < -pi] <- dxdy[dxdy < -pi] + 2*pi
   dxdy <- dxdy/2
   dxdy[c(1, N)] <- NA
   dxdy
}


ReIm <- function(complex) {
   list(R = Re(complex), I = Im(complex))
}

clusters <- function(data, k) {
   kmeans(data, k)$cluster
}

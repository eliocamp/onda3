world <- subset(map_data("world2"), lat > 0)

geom_world <- geom_path(data = world, aes(long, lat, group = group),
                        size = 0.2)



formula <- u ~ year | lon + lat

Regression2D <- function(formula, data = NULL) {
  f <- as.character(formula)
  f <- stringr::str_split(f,"~", n = 2)[[1]]
  
  value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])
  
  matrix.vars <- f[stringr::str_detect(f, "\\|")]
  matrix.vars <- stringr::str_split(matrix.vars,"\\|", n = 2)[[1]]
  
  row.vars <- stringr::str_squish(stringr::str_split(matrix.vars[1], "\\+")[[1]])
  col.vars <- stringr::str_squish(stringr::str_split(matrix.vars[2], "\\+")[[1]])
  
  if (is.null(data)) {
    formula <- Formula::as.Formula(formula)
    data <- as.data.table(eval(quote(model.frame(formula, data  = data))))
  } else {
    # Check if columns are indata
    all.cols <- c(value.var, row.vars, col.vars)
    missing.cols <- all.cols[!(all.cols %in% colnames(data))]
    if (length(missing.cols) != 0) {
      stop(paste0("Columns not found in data: ", paste0(missing.cols, collapse = ", ")))
    }
    data <- setDT(data)[, (all.cols), with = FALSE]
  }
  
  setDT(data)
  dcast.formula <- stringr::str_squish(f[stringr::str_detect(f, "\\|")])
  dcast.formula <- as.formula(stringr::str_replace(dcast.formula, "\\|", "~"))
  value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])
  
  g <- metR:::.tidy2matrix(data, dcast.formula, value.var, fill = NULL)
  
  if (length(g$matrix) < nrow(data)) {
    stop(paste("The formula ", as.character(formula), " does not identify an unique observation for each cell."))
  }
  
  N <- nrow(g$matrix)
  eof <- svd(g$matrix)
  set.seed(42)
  fit <- glmnet::cv.glmnet(eof$u*sqrt(N), g$rowdims[[1]])
  
  p <- predict(fit$glmnet.fit, newx = eof$u, s = fit$lambda.1se)
  
  coef <- c(as.matrix(coef(fit, s = fit$lambda.1se)))[-1]
  
  E <- eof$v%*%diag(eof$d, nrow = length(eof$d))/sqrt(N)
  coef <- t(matrix(coef)) %*% t(E)
  
  set(g$coldims, NULL, value.var, c(coef))
  r2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == fit$lambda.1se)] 
  
  set(g$coldims, NULL, "r2", c(r2))
  return(g$coldims)
}


data <- ReadNetCDF("~/DATOS/NCEP Reanalysis/uwnd.mon.mean.nc", c(u = "uwnd"),
                   subset = list(level = 300,
                                 lat = c(0:90),
                                 time = c("1948-01-01", "2009-12-31"))) %>% 
  .[month(time) %in% c(12, 1, 2)] %>% 
  .[, .(u = mean(u)), by = .(lat, lon, year(time))] 

data[, u := Anomaly(u), by = .(lat, lon)]


resample <- function(data) {
  years <- sample(unique(data$year), replace = FALSE)
  return(copy(data)[, year := years, by = .(lat, lon)])
}

point_regression <- data[, FitLm(u, year), by = .(lon, lat)]

ggplot(point_regression[term == "year"], aes(lon, lat)) +
  geom_contour_fill(aes(z = estimate*10)) +
  geom_world +
  scale_fill_divergent() +
  coord_quickmap()



data[, Regression2D(u ~ year | lat + lon)] %>% 
  ggplot(aes(lon, lat)) +
  geom_contour_fill(aes(z = u*10)) +
  geom_world +
  scale_fill_divergent() +
  coord_quickmap()

B <- 100
fdr_naive <- unlist(lapply(1:B, function(x) {
  resample(data) %>% 
    .[, FitLm(u, year, se = TRUE), by = .(lon, lat)] %>% 
    .[term == "year"] %>%
    .[, mean(abs(estimate)/std.error > 2)]
}))

data <- resample(data)

data %>% 
  .[, FitLm(u, year, se = TRUE), by = .(lon, lat)] %>% 
  .[term == "year"] %>%
  ggplot(aes(lon, lat)) +
  geom_contour_fill(aes(z = estimate*10)) +
  geom_world +
  scale_fill_divergent() +
  coord_quickmap()

data %>% 
  .[, Regression2D(u ~ year | lon + lat)] %>% 
  ggplot(aes(lon, lat)) +
  geom_contour_fill(aes(z = u*10)) +
  geom_world +
  scale_fill_divergent() +
  coord_quickmap()



r2_eof <- unlist(lapply(1:B, function(x) {
  resample(data) %>% 
    Regression2D(u ~ year | lon + lat, data = .) 
}))


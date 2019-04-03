library(ggplot2)
library(metR)
library(data.table)
library(magrittr)

world <- subset(map_data("world2"), lat > 0)

geom_world <- geom_path(data = world, aes(long, lat, group = group),
                        size = 0.2)

Regression2D <- function(formula, data = NULL, B = 1, seed = 42) {
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
   N <- length(g$rowdims[[1]])
   
   
   eof <- svd(g$matrix)
   fit <- withr::with_seed(seed, glmnet::cv.glmnet(eof$u*sqrt(N), g$rowdims[[1]], standardize = FALSE))
   lambda.min <- fit$lambda.min
   r2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == fit$lambda.min)]
   coef <- c(as.matrix(coef(fit, s = lambda.min)))[-1]
   
   nonzero <- coef != 0
   fit <- .lm.fit(cbind(1, eof$u[, nonzero]*sqrt(N)), g$rowdims[[1]])
   coef <- coef(fit)[-1]
   
   r2 <- 1 - var(fit$residuals)/var(g$rowdims[[1]])
   # k <- 21
   # eof <- svd(g$matrix, k, k)
   # eof$d <- eof$d[1:k] 
   # coef <- coef(.lm.fit(cbind(1, eof$u*sqrt(N)), g$rowdims[[1]]))[-1]
   
   E <- eof$v%*%diag(eof$d, nrow = length(eof$d))/sqrt(N)
   coef <-  coef %*% t(E[, nonzero]) / var(g$rowdims[[1]])
   
   set(g$coldims, NULL, value.var, c(coef))
   
   if (B > 1) {
      r2_null <- withr::with_seed(seed, {
         vapply(seq_len(B), function(b) {
            Prow <- sample(seq_len(N), replace = TRUE)
            eof <- svd(g$matrix[Prow, ])
            fit <- glmnet::cv.glmnet(eof$u, g$rowdims[[1]], standardize = FALSE)
            fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == fit$lambda.min)] 
         }, 1)})
      r2_c <- mean(r2_null >= r2)
   } else {
      r2_c <- NULL
   }
   
   # setattr(g$coldims, "r2", r2)
   # setattr(g$coldims, "p.value", r2_c)
   # set(g$coldims, NULL, "r2", c(r2))
   # set(g$coldims, NULL, "r2_c", c(r2_c))
   
   # return(list(field = g$coldims,
   #             r2 = r2,
   #             p.value = r2_c))
   g$coldims
}


data <- ReadNetCDF("~/DATOS/NCEP Reanalysis/uwnd.mon.mean.nc", c(u = "uwnd"),
                   subset = list(level = 300,
                                 lat = c(0:90),
                                 time = c("1948-01-01", "2009-12-31"))) %>% 
   .[month(time) %in% c(12, 1, 2)] %>% 
   .[, .(u = mean(u)), by = .(lat, lon, year(time))] 

data[, u := Anomaly(u), by = .(lat, lon)]

y <- unique(data$year)
data <- data[year %in% y[1:n]]

resample <- function(data) {
   years <- sample(unique(data$year), replace = FALSE)
   return(copy(data)[, year := years, by = .(lat, lon)])
}

point_regression <- data[, FitLm(u, year, se = TRUE), by = .(lon, lat)]

ggplot(point_regression[term == "year"], aes(lon, lat)) +
   geom_contour_fill(aes(z = estimate)) +
   geom_world +
   scale_fill_divergent("Trend in U(300hpa)") +
   coord_quickmap() +
   metR:::theme_field()



data[, Regression2D(u ~ year | lat + lon)] %>% 
   ggplot(aes(lon, lat)) +
   geom_contour_fill(aes(z = u)) +
   geom_world +
   scale_fill_divergent("Trend in U(300hpa)") +
   coord_quickmap() +
   metR:::theme_field()


data[, FitLm(u, year, se = TRUE), by = .(lon, lat)] %>% 
   .[term == "year"] %>% 
   .[data[, Regression2D(u ~ year | lat + lon)], on = c("lon", "lat")] %>% 
   ggplot(aes(u, estimate)) +
   geom_point() +
   geom_abline() +
   coord_fixed()



y <- unique(data$year)
new_data <- data[year %in% y[1:n]]
ratio <- unlist(lapply(10:N, function(n) {
   new_data <- data[year %in% y[1:n]]
   new_data[, FitLm(u, year, se = TRUE), by = .(lon, lat)] %>% 
      .[term == "year"] %>% 
      .[new_data[, Regression2D(u ~ year | lat + lon)$field], on = c("lon", "lat")] %>% 
      lm(u ~ estimate, data = .) %>% 
      coef() %>% 
      .[["estimate"]]   
}

))

ggplot(data.frame(N = 10:N, ratio = ratio), aes(N^3, ratio)) +
   geom_point() +
   geom_smooth(method = "lm")

coef(lm(ratio ~ I(N^3) - 1, data = data.frame(N = 10:N, ratio = ratio)))


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
      Regression2D(u ~ year | lon + lat, data = ., B = 100) %>% 
      attr("p.value")
}))


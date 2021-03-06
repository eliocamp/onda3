Regression2D <- function(formula, y, data = NULL, 
                         method = c("cv", "lasso", "neof"),
                         max_eof = Inf,
                         k_fold = 10,
                         alpha = 1, 
                         seed = 42) {
   method <- method[1]
   if (method == "neof") {
      fit_function <- function(X, y) {
         zeros <- rep(0, length = ncol(X) - max_eof)
         
         X <- X[, seq_len(max_eof), drop = FALSE]
         fit <- .lm.fit(cbind(1, X), y)
         r2 <- 1 - var(fit$residuals)/var(y)
         
         return(list(coef = c(coef(fit)[-1], zeros),
                     r2 = r2))
      }
   } else if (method == "lasso") {
      fit_function <- function(X, y) {
         fit <- withr::with_seed(seed, glmnet::cv.glmnet(X, y, nfolds = k_fold,
                                                         lambda.min.ratio = 0.01,
                                                         alpha = alpha,
                                                         standardize = FALSE))
         lambda <- fit$lambda.1se
         r2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == lambda)]
         coef <- c(as.matrix(coef(fit, s = lambda)))[-1]
         return(list(r2 = r2, coef = coef))
      }
   } else if (method == "cv") {
      
      fit_function <- function(X, y) {
         n_features <- ncol(X)
         
         if (k_fold > N) {
            warning("K-folding for crossvalidation greater than N.", 
                    " Reducing to leave-one-out crossvalidation")
            k_fold <- N
         }
         
         N <- length(y)
         start <- floor(seq(1, N, by = N/k_fold))
         end   <-  c((start - 1)[-1], N)
         
         mean_error <- rep(0, length = n_features)
         sd_error <-  rep(0, length = n_features)
         for (f in seq_len(n_features)) {
            error <- rep(0, N)
            for (k in seq_along(k_fold)) {
               
               test <- seq.int(start[k], end[k])
               
               fit <- .lm.fit(cbind(1, X[-test, seq_len(f), drop = FALSE]), y[-test])
               error[test] <-  cbind(1, X[test, seq_len(f), drop = FALSE]) %*% coef(fit) - y[test]
            }
            mean_error[f] <- mean(error^2)
            sd_error[f] <- sd(error^2)/sqrt(length(error))
         }
         
         min_error <- which.min(mean_error)
         
         features <- which((mean_error - sd_error)[seq_len(min_error)] < (mean_error + sd_error)[min_error])[1]
         fit <- .lm.fit(cbind(1, X[, seq_len(features), drop = FALSE]), y)
         coef <- c(coef(fit)[-1], rep(0, length = ncol(X) - features))
         r2 <- 1 - var(fit$residuals)/var(y)
         
         return(list(coef = coef,
                     r2 = r2))
      }
   } else {
      stop("unknown method")
   }
   
   
   # Housekeeping. Getting data and transforming it to matrix
   formula <- enrich_formula(formula)
   data <- data_from_formula(formula = formula, 
                             data = data)
   g <- metR:::.tidy2matrix(data, formula$dims, formula$value.var, fill = NULL)
   
   if (length(g$matrix) < nrow(data)) {
      stop(paste("The formula", as.character(formula), "does not identify an unique observation for each cell."))
   }
   
   N <- length(g$rowdims[[1]])
   
   force(y)
   y_name <- deparse(substitute(y))
   y <- metR:::.tidy2matrix(data, formula$dims, y_name, fill = NULL)$matrix[, 1]
   
   g$matrix <- scale(g$matrix, scale = FALSE)
   
   if (max_eof < 1) {
      max_eof <- round(max_eof*min(nrow(g$matrix), ncol(g$matrix)))
   } else {
      max_eof <- min(max_eof, nrow(g$matrix) - 1, ncol(g$matrix) - 1) 
   }
   
   if (max_eof >= 0.5*min(nrow(g$matrix), ncol(g$matrix))) {
      svd_fun <- base::svd
   } else {
      svd_fun <- irlba::irlba
   }
   
   eof <- svd_fun(g$matrix, max_eof, max_eof)
   eof$d <- eof$d[seq_len(max_eof)]
   
   
   fit <- fit_function(eof$u*sqrt(N), y)
   non_zero <- which(fit$coef != 0)
   
   E <- eof$v%*%diag(eof$d, nrow = length(eof$d))/sqrt(N)
   coef_eof <- fit$coef
   fit$coef <-  fit$coef %*% t(E) / var(y)   # Why divide by variance of y?
   
   M <- length(non_zero)
   N <- length(y)
   f.statistic <-  fit$r2/(1 - fit$r2)*(N - M - 1)/M
   set(g$coldims, NULL, formula$value.var, c(fit$coef))
   # set(g$coldims, NULL, "p.value", c(p.val_field))
   return(
      list(field = g$coldims,
           coef_eof = coef_eof,
           summary = data.frame(r2 = fit$r2,
                                f.statistic = f.statistic,
                                p.value = pf(f.statistic, N, N - M - 1, lower.tail = FALSE),
                                non_zero = sum(coef_eof != 0)
           )
      ))
}



enrich_formula <- function(formula) {
   f <- as.character(formula)
   f <- stringr::str_split(f,"~", n = 2)[[1]]
   
   value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])
   
   matrix.vars <- f[stringr::str_detect(f, "\\|")]
   matrix.vars <- stringr::str_split(matrix.vars,"\\|", n = 2)[[1]]
   
   row.vars <- stringr::str_squish(stringr::str_split(matrix.vars[1], "\\+")[[1]])
   col.vars <- stringr::str_squish(stringr::str_split(matrix.vars[2], "\\+")[[1]])
   
   dcast.formula <- stringr::str_squish(f[stringr::str_detect(f, "\\|")])
   dcast.formula <- as.formula(stringr::str_replace(dcast.formula, "\\|", "~"))
   value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])
   
   return(list(formula = formula,
               dims = dcast.formula, 
               value.var = value.var,
               row.vars = row.vars, 
               col.vars = col.vars))
}

data_from_formula <- function(formula, data) {
   if (is.null(data)) {
      data <- as.data.table(eval(quote(model.frame(Formula::as.Formula(formula$formula),
                                                   data  = data))))
   } else {
      # Check if columns are indata
      all.cols <- c(formula$value.var, formula$row.vars, formula$col.vars)
      missing.cols <- all.cols[!(all.cols %in% colnames(data))]
      if (length(missing.cols) != 0) {
         stop(paste0("Columns not found in data: ", paste0(missing.cols, collapse = ", ")))
      }
      data <- setDT(data)[, (all.cols), with = FALSE]
   }
   return(data)
}

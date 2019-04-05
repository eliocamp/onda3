
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

lasso_fit <- function(X, y, seed) {
   fit <- withr::with_seed(seed, glmnet::cv.glmnet(X, y, standardize = FALSE))
   lambda <- fit$lambda.1se
   r2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == lambda)]
   coef <- c(as.matrix(coef(fit, s = lambda)))[-1]
   return(list(r2 = r2, coef = coef))
}

Regression2D <- function(formula, y, data = NULL, B = 1, 
                         method = c("cv", "lasso", "neof"),
                         keep = 0.95, seed = 42) {
   # Housekeeping. Getting data and transforming it to matrix
   formula <- enrich_formula(formula)
   data <- data_from_formula(formula = formula, 
                             data = data)
   g <- metR:::.tidy2matrix(data, formula$dims, formula$value.var, fill = NULL)
   
   if (length(g$matrix) < nrow(data)) {
      stop(paste("The formula ", as.character(formula), " does not identify an unique observation for each cell."))
   }
   
   N <- length(g$rowdims[[1]])
   if (missing(y)) y <- as.numeric(g$rowdims[[1]])
   
   # Algo: get EOF, fit on principal components, multiply with fields to 
   # get regression field. 
   
   g$matrix <- scale(g$matrix, scale = FALSE)
   eof <- svd(g$matrix)
   
   if (method == "n_eof") {
      # TEST
      # Limit to first EOF that explain 95% variance en be done with it. 
      v.g <- norm(abs(g$matrix), type = "F")
      exp_var <- eof$d^2/v.g^2
      keep <- which(cumsum(exp_var) >= keep)[1]
      
      fit_function <- function(X, y, seed) {
         zeros <- rep(0, length = ncol(X) - keep)
         X <- X[, seq_len(keep), drop = FALSE]
         fit <- .lm.fit(cbind(1, X), y)
         r2 <- 1 - var(fit$residuals)/var(y)
         
         return(list(coef = c(coef(fit)[-1], zeros),
                     r2 = r2))
      }
   } else if (method == "lasso") {
      fit_function <- function(X, y, seed) {
         fit <- withr::with_seed(seed, glmnet::cv.glmnet(X, y, standardize = FALSE))
         lambda <- fit$lambda.1se
         r2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == lambda)]
         coef <- c(as.matrix(coef(fit, s = lambda)))[-1]
         return(list(r2 = r2, coef = coef))
      }
   } else if (method == "cv") {
      fit_function <- function(X, y, seed) {
         K_max <- min(ncol(X), length(y))
         mean_error <- rep(0, length = K_max)
         sd_error <-  rep(0, length = K_max)
         for (K in seq_len(K_max)) {
            error <- rep(0, length = length(y))
            for (n in seq_along(y)) {
               fit <- .lm.fit(cbind(1, X[-n, seq_len(K), drop = FALSE]), y[-n])
               error[n] <-  cbind(1, X[n, seq_len(K), drop = FALSE]) %*% coef(fit) - y[n]
            }
            mean_error[K] <- mean(error^2)
            sd_error[K] <- sd(error^2)/sqrt(length(error))
         }
         
         min_error <- which.min(mean_error)
         
         K <- which(mean_error[seq_len(min_error)] < (mean_error + sd_error)[min_error])[1]
         fit <- .lm.fit(cbind(1, X[, seq_len(K), drop = FALSE]), y)
         coef <- c(coef(fit)[-1], rep(0, length = ncol(X) - K))
         r2 <- 1 - var(fit$residuals)/var(y)
         
         return(list(coef = coef,
                     r2 = r2))
      }
   } else {
      stop("unknown method")
   }

   fit <- fit_function(eof$u*sqrt(N), y, seed)
   non_zero <- which(fit$coef != 0)
   
   E <- eof$v%*%diag(eof$d, nrow = length(eof$d))/sqrt(N)
   fit$coef <-  fit$coef %*% t(E) / var(y)   # Why divide by variance of y?
   
   # Bootstrap
   # if (B > 1 & any(fit$coef != 0)) {
   #    withr::with_seed(seed, {
   #       p.val_field <- vector("numeric", length = length(fit$coef))
   #       p.val_r2 <- 0
   #       for (b in seq_len(B)) {
   #          scramble <- sample(seq_len(N), replace = TRUE)
   #          
   #          new_fit <- fit_function(eof$u*sqrt(N), y[scramble], seed)
   #          new_fit$coef <-  new_fit$coef %*% t(E) / var(y) 
   # 
   #          
   #          new_p <- (abs(new_fit$coef) >= abs(fit$coef))
   #          p.val_field <- p.val_field + new_p
   #          
   #          p.val_r2    <- p.val_r2 + (new_fit$r2 >= fit$r2)
   #       }
   #       
   #       p.val_r2 <-  p.val_r2/B
   #       p.val_field <-  p.val_field/B
   #    })
   # } else {
   #    p.val_r2 <-  NA
   #    p.val_field <- NA
   # }
   
   M <- length(non_zero)
   N <- length(y)
   f.statistic <-  fit$r2/(1 - fit$r2)*(N - M - 1)/M
   set(g$coldims, NULL, formula$value.var, c(fit$coef))
   # set(g$coldims, NULL, "p.value", c(p.val_field))
   return(
      list(field = g$coldims,
               non.zero = non_zero,
               r2 = fit$r2,
               f.statistic = f.statistic,
               p.value = pf(f.statistic, N, N - M - 1, lower.tail = FALSE) 
   ))
}


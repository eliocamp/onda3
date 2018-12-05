

sim_stationarity <- function(n, d = 0, trans = FALSE) {
   x <- rnorm(n) + d
   y <- rnorm(n)
   am <- sqrt(x^2 + y^2)
   phase <- atan2(y, x)
   if (isTRUE(trans)) {
      2/pi*asin(weighted.mean(cos(phase - mean.phase(am, phase, 1)), am))
   } else {
      weighted.mean(cos(phase - mean.phase(am, phase, 1)), am)
   }
}

modelise <- function(S, n) {
   model <- nls(S ~ b*n^g + S0, start = c(S0 = 0.1, g = -1/2, b = 1))
   as.list(coef(model))
}


B <- 500
N <- seq(2, 300, by = 20)
ds <- 10^seq(-2, 0.5, length.out = 20)

trials <- data.table(expand.grid(n = N,
                                 d = ds,
                                 trans = c(TRUE, FALSE)))

trials <- trials[, .(S = mean(sapply(1:B, function(x) sim_stationarity(n, d, trans)))),
                 by = .(n, d, trans)]


coefs <- trials[, modelise(S, n), by = .(d, trans)]

pred <- coefs[,  .(n = N, 
           S = b*N^g + S0), by = .(d, trans)] 

ggplot(trials, aes(n, S, color = d, group = d)) +
   geom_point() +
   geom_line(data = pred) +
   facet_wrap(~trans) 

ggplot(coefs, aes(S0, b)) +
   geom_point() +
   stat_function(fun = function(x) (1 - x)) +
   facet_wrap(~trans)

ggplot(coefs, aes(S0, g)) +
   geom_point() +
   stat_function(fun = function(x) (-1/2 - x)) +
   facet_wrap(~trans)





x <- rnorm(30)
y <- rnorm(30)
am <- sqrt(x^2 + y^2)
phase <- atan2(y, x)

.S <- function(am, phase, n = length(x)){
   select <- sample(seq_along(am), n, replace = FALSE)
   phase <- phase[select]
   am <- am[select]
   
   2/pi*asin(weighted.mean(cos(phase - mean.phase(am, phase, 1)), am))
}

.Ss <- function(am, phase) {
   s_hat <- sapply(seq_along(am)[-1], function(i) .S(am, phase, n = i))
   n <- seq_along(am)[-1]
   # plot(n, s_hat)
   model <- suppressWarnings(nls(s_hat ~ (1-S0)*n^g + S0, start = c(g = -1/2, S0 = 0.1), 
                control = nls.control(warnOnly = TRUE)))
   if(model$convInfo$isConv == FALSE) {
      return(NA)
   }
   
   unname(coef(model)[2])
}

Ss <- function(am, phase, B = 10) {
   s_hat <- c(sapply(seq_len(B), function(b) .Ss(am, phase)))
   s_hat <- s_hat[!is.na(s_hat)]
   S <- mean(s_hat)
   std.err <- sd(s_hat)/sqrt(length(s_hat))
   S - std.err*2 < 0 
}


.Ss2 <- function(am, phase) {
   sapply(seq_along(am)[-1], function(i) .S(am, phase, n = i))
}

Ss2 <- function(am, phase, B = 10) {
   s_hat <- c(sapply(seq_len(B), function(b) .Ss2(am, phase)))
   n <- rep(seq_along(am)[-1], B)
   model <- nls(s_hat ~ (1-S0)*n^g + S0, start = c(g = -1/2, S0 = 0.1))
   a <- summary(model)$coefficients[2, 1:2]
   unname(a[1] - a[2]*2 < 0)
}


# tests


x <- rnorm(30)
y <- rnorm(30)
am <- sqrt(x^2 + y^2)
phase <- atan2(y, x)

B <- 200
tests <- 100

tests <- rbindlist(lapply(1:tests, 
                         function(t) {
                            x <- rnorm(30)
                            y <- rnorm(30)
                            am <- sqrt(x^2 + y^2)
                            phase <- atan2(y, x)
                            
                            one_true_model <- Ss2(am, phase, B = B)
                            small_models <- Ss(am, phase, B = B)
                            data.table(one_true_model, small_models)
                         }))

melt(tests) %>% 
   .[, mean(value), by = variable]



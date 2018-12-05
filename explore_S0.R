

stat1 <- function(n, d = 1) {
   x <- rnorm(n) + d
   y <- rnorm(n)
   am <- sqrt(x^2 + y^2)
   phase <- atan2(y, x)
   
   weighted.mean(cos(phase - mean.phase(am, phase, 1)), am)
}

modelise <- function(S, n) {
   model <- nls(S ~ (1-S0)*n^g + S0, start = c(g = -1/2, S0 = 0.5))
   as.list(coef(model))
}


B <- 400
N <- seq(1, 1580, length.out = 60)
ds <- 0

sims <- as.data.table(expand.grid(n = N))

sims <- sims[, .(k = sapply(1:B, function(x) stat1(n, d = 0))), by = .(n)]

sims[, S := 2/pi*asin(k)]

S0 <- melt(sims, id.vars = c("n")) %>% 
   .[, modelise(value, n), by = variable]

melt(sims, id.vars = c("n")) %>% 
   .[, mean(value), by = .(n, variable)] %>% 
   ggplot(aes(n, V1)) +
   geom_line(aes(color = variable))

N <-N <- seq(1, 500000, length.out = 60)
N^(-1/2)/(N^(-0.68) + 0.007)

ggplot(sims[d == 0], aes(n, S)) +
   geom_line() +
   stat_function(fun = function(x) x^(-1/2))

ggplot(S0, aes(S0, g)) +
   geom_point() 
   stat_function(fun = function(x) 1 - x)

ggplot(sims, aes(n, S)) +
   geom_line(aes(color = d, group = d)) +
   geom_point() 

coefs <- sims[, modelise(S, n), by = d]

   ggplot(coefs, aes(d, S0)) +
   geom_point()  +
   stat_function(fun = function(x) 1.963/(1 + exp(-1.65*x))-1)
   # geom_smooth(method = "lm")

model <-  nls(S0 ~ L/(1 + exp(k*d)) - 1, data = coefs, start = c(L = 1, k = 1))




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



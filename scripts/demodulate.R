period <- 20 #360/3
fc <- 1/period
fs <- 1
n <- 360
t <- seq(0, 360)
AMP <- 1.2 + sin(t*pi/180)
PH  <- 5 * t / max(t)
noise <- 0 # rnorm(n, sd=0.5)
carrier <- sin(2 * pi * fc * t + PH*pi/180)
signal <- AMP * carrier
x <- signal + noise
# 
x <- vwnd[time == time[1] & hemisphere == "SH", v]
t <- unique(vwnd$lon)
fs <- 1/ggplot2::resolution(unique(vwnd$lon))

xc <- x * cos(2 * pi * fc * t)
xs <- x * sin(2 * pi * fc * t)
par(mar=c(1.75,3,1/2,1), mgp=c(2, 0.7, 0), mfcol=c(3,1))
plot(t, x, type='l')
plot(t, cos(2 * pi * fc * t), type='l')
plot(t, xc, type='l')


w <- 2 * fc / fs
## Here, we use more smoothing
w <- w * 1.
filter <- signal::butter(5, w)    # FIXME: why extras on w?
xcs <- signal::filtfilt(filter, xc)
xss <- signal::filtfilt(filter, xs)
amp <- 2 * sqrt(xcs^2 + xss^2)
phase <- 180 / pi * atan2(xcs, xss)

par(mar=c(1.75,3,1/2,1), mgp=c(2, 0.7, 0), mfcol=c(3,1))
plot(t, amp, type='l')
lines(t, amp, col='red')

plot(t, phase, type='l')
lines(t, phase, col='red')
plot(t, x, type='l', pch=20)
lines(t, amp * sin(2 * pi * fc * t + phase*pi/180), col='red')


Demodulate <- function(x, t, k = 3, mult = seq(0.1, 5, length.out = 10)) {
   # x <- vwnd[time == time[1] & hemisphere == "SH", v]
   # t <- unique(vwnd$lon)
   x_m <- mean(x)
   x_a <- x - x_m
   
   fs <- 1/ggplot2::resolution(t, zero = FALSE)
   
   var_x <- var(x_a)
   # Avoid edge effects
   keep <- (length(x_a)+1):(2*length(x_a))
   x_a <- rep(x_a, 3)
 
   # from https://dankelley.github.io/r/2014/02/17/demodulation.html
   ks <- k
   demod <- lapply(ks, function(k) {
      fc <- k/360
      t <- c(t-360, t, t+360)
      xc <- x_a * cos(k*pi/180 * t)
      xs <- x_a * sin(k*pi/180 * t)
      w_base <- 2 * fc / fs
      if (length(mult) > 1) {
         r2s <- lapply(mult, function(m) {
            ## Here, we use more smoothing
            w <- w_base 
            filter <- signal::butter(5, w)    # FIXME: why extras on w?
            xcs <- signal::filtfilt(filter, xc)
            xss <- signal::filtfilt(filter, xs)
            amplitude <- 2 * sqrt(xcs^2 + xss^2)
            phase <- -(180 / pi * atan2(xcs, xss) - 180/2) /k
            carrier <- cos(k*pi/180 * (t - phase))
            pred <- amplitude * carrier
            
            r2 <- 1 - var((x_a-pred)[keep])/var_x
            
            return(list(r2 = r2, 
                        amplitude = amplitude[keep], 
                        carrier = carrier[keep]))
         })
         
         max <- which.max(vapply(r2s, function(x) x$r2, 1))
         
         return(with(r2s[[max]], 
                     list(amplitude = amplitude,
                          carrier = carrier,
                          offset = x_m)))
      } else {
         # browser()
         w <- w_base 
         filter <- signal::butter(5, w)    # FIXME: why extras on w?
         xcs <- signal::filtfilt(filter, xc)
         xss <- signal::filtfilt(filter, xs)
         amplitude <- 2 * sqrt(xcs^2 + xss^2)
         phase <- -(180 / pi * atan2(xcs, xss) - 180/2) /k
         carrier <- cos(k*pi/180 * (t - phase))
         
         return(list(amplitude = amplitude[keep], 
                     carrier = carrier[keep],
                     offset = x_m))
      }
   })
   # browser()
   return(demod)
   
}


demodulate_one <- function(xc, xs, t, k, w) {

   return(list(amplitude = amplitude, 
               phase = phase))
}

vwnd[, c("ampl3", "phase3") := Demodulate(v, lon, 3, 1), by = .(time, hemisphere)]







t <- unique(vwnd$lon)

amp1 <- 1.55 + cos(t*pi/180 - 30*pi/180) + 0.2*cos(2*t*pi/180 - 15*pi/180)
amp <- amp1 +  0.6*cos(5*t*pi/180)
carrier <- 1*cos(3*(t - 10)*pi/180)
x <- carrier*amp 



k <- 3



Demodulate  <- function(x, k) {
   t <- seq(0, 360 - 360/length(x), by = 360/length(x))*pi/180
   
   # Onda carrier
   fit1 <- FitWave(x, k = k)
   carrier <- cos(k*(t - fit1$phase))
   
   # "bases" de senos y cosenos de la amplitud
   js <- seq_len(k-1)
   cosjs <- lapply(js, function(j) cos(t*j)*carrier)
   sinjs <- lapply(js, function(j) sin(t*j)*carrier)
   
   # Fit de todo esto
   X <- do.call(cbind, c(list(carrier), cosjs, sinjs))
   new_fit <- .lm.fit(X, FilterWave(x, 0:3))
   
   A <- (X/carrier) %*% new_fit$coefficients
   
   return(list(carrier = carrier, 
               amplitude = A))
}
```





plot(amp1)
env <- Demodulate(x, 3)
lines((env - mean(env))*mean(env) + mean(env))

plot(x)
lines(FilterWave(x, 1:3))

plot(amp1)
env <- WaveEnvelope(FilterWave(x, 0:3))
lines(env)


plot(amp1)
lines(env)

x <- vwnd[time == unique(time)[20] & hemisphere == "SH", v]
t <- unique(vwnd$lon)
A <- Demodulate(x, t, 3)
zw3 <- FilterWave(x, 3)
plot(x)
lines(A*zw3)
lines(zw3)
lines(FilterWave(x, 0:3), col = "red")
lines(FilterWave(x, 0:3, -1))

Demodulate(x, 3)

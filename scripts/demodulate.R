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
   x <- vwnd[time == time[1] & hemisphere == "SH", v]
   t <- unique(vwnd$lon)
   
   fs <- 1/ggplot2::resolution(t, zero = FALSE)
   
   var_x <- var(x)
   # Avoid edge effects
   keep <- (length(x)+1):(2*length(x))
   x <- rep(x, 3)
   t <- c(t-360, t, t+360)

   # from https://dankelley.github.io/r/2014/02/17/demodulation.html

   ks <- k
   demod <- lapply(ks, function(k) {
      fc <- k/360
      
      xc <- x * cos(k*pi/180 * t)
      xs <- x * sin(k*pi/180 * t)
      w_base <- 2 * fc / fs
      if (length(mult) > 1) {
         r2s <- lapply(mult, function(m) {
            ## Here, we use more smoothing
            w <- w_base * m
            filter <- signal::butter(5, w)    # FIXME: why extras on w?
            xcs <- signal::filtfilt(filter, xc)
            xss <- signal::filtfilt(filter, xs)
            amplitude <- 2 * sqrt(xcs^2 + xss^2)
            phase <- -(180 / pi * atan2(xcs, xss) - 180/2) /k
            
            pred <- amplitude * cos(k*pi/180 * (t - phase))
            
            r2 <- 1 - var((x-pred)[keep])/var_x
            
            return(list(r2 = r2, 
                        amplitude = amplitude[keep], 
                        phase = phase[keep]))
         })
         
         max <- which.max(vapply(r2s, function(x) x$r2, 1))
         
         return(with(r2s[[max]], 
                     list(amplitude = amp,
                          phase = phase)))
      } else {
         w <- w * mult
         filter <- signal::butter(5, w)    # FIXME: why extras on w?
         xcs <- signal::filtfilt(filter, xc)
         xss <- signal::filtfilt(filter, xs)
         amp <- 2 * sqrt(xcs^2 + xss^2)
         phase <- 180 / pi * atan2(xcs, xss)
         return(list(amplitude = amp[keep], 
                     phase = phase[keep]))
      }
   })
   
}


demodulate_one <- function(xc, xs, t, k, w) {

   return(list(amplitude = amplitude, 
               phase = phase))
}

vwnd[, c("ampl3", "phase3") := Demodulate(v, lon, 3, 1), by = .(time, hemisphere)]





Demodulate  <- function(x, k) {
   fit1 <- FitWave(x, k = k)
   
   cosfit <- cos(k*(t*pi/180 - fit1$phase))
   
   A <- FilterWave(x, 1:3)/cosfit
   
   near_zero <- abs(cosfit) < 0.0001
   
   A_full <- approx(seq_along(x)[!near_zero], A[!near_zero], xout = seq_along(x))
   
   FilterWave(A_full$y, seq(0, k-1))
}



t <- unique(vwnd$lon)

amp1 <- 1.5 + cos(t*pi/180 - 30*pi/180) 
amp <- amp1 + 0.3*cos(5*t*pi/180)
carrier <- 1*cos(3*(t - 10)*pi/180)
x <- carrier*amp 


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

x <- vwnd[time == time[1] & hemisphere == "SH", v]

Demodulate(x, 3)

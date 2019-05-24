res <- 2.5
t <-  seq(0, 360 - res, by = res)

# Modelación de la amplitud (onda 1 y onda 2)
amp <- 1.5 + 0.5*cos(4*t*pi/180 - 30*pi/180) + rnorm(length(x)) #+ 0.6*cos(2*t*pi/180 - 25*pi/180) 


# Onda 3 pura
carrier <- cos(3*(t - 10)*pi/180) 

# Señal total 
x <- carrier*amp   #+ 0.8*cos(5*t*pi/180)

f <- signal::specgram(x, n = length(x)/3, Fs = res, overlap = length(x)/3 - 1)

f$f <- length(x)*(f$f)/res

by =  diff(range(t))/(length(f$t))
ts <- seq(t[1] + by/2, t[length(t)] - by/2, by = by)
f$t <- ts
dimnames(f$S) <- list(f = f$f, t = f$t)



f.df <- as.data.table(melt(abs(f$S)/nrow(f$S)*2))


ggplot(f.df[f == 3], aes(t, value)) +
   geom_line(data = data.frame(t = t, value = x), color = "gray") +
   geom_line() +
   geom_point()  +
   geom_line(data = data.frame(t = t, value = amp), color = "red") +
   # geom_line(data = data.frame(t = t, value = WaveEnvelope(x))) +
   geom_line(data = data.frame(t = t, value = WaveEnvelope(x))) +
   # geom_line(data = data.frame(t = t, value = PeriodicWavelet(x, 3)[[1]]))  +
   NULL


fortify.list <- function(model, data, ... ){
   as.data.frame(model)
}

x <- ej$hgt.z

as.data.frame(FitWave(x, 1:10)) %>% 
   ggplot(aes(k, amplitude)) +
   geom_line()


ggplot(f.df[f == 3], aes(t, value)) +
   geom_line(data = data.frame(t = t, value = x), color = "gray") +
   # geom_line() +
   # geom_point()  +
   # geom_line(data = data.frame(t = t, value = amp), color = "red") +
   # geom_line(data = data.frame(t = t, value = WaveEnvelope(x))) +
   geom_line(data = data.frame(t = t, value = WaveEnvelope(x, 1:3))) +
   geom_line(data = data.frame(t = t, value = WaveEnvelope(x, 1:3)*FilterWave(x, 2)/FitWave(x, 2)$amplitude))
# geom_line(data = data.frame(t = t, value = PeriodicWavelet(x, 3)[[1]]))  +
NULL




SH %>% 
   .[, FitWave(hgt, 1:6), by = time] %>% 
   .[, .(k = k[which.max(amplitude)]), by =.(time)] -> main_wave

SH <- main_wave[SH, on = "time"]

SH[,  hgt_filt := FilterWave(hgt, seq(1, k[1])), by = time]
SH[, hgt_base := FilterWave(hgt, k[1])/FitWave(hgt, k[1])$amplitude, by = time]

SH[, env := WaveEnvelope(hgt_filt), by = time]



SH[time == unique(time)[33]] %>% 
   ggplot(aes(lon, Anomaly(hgt))) +
   geom_line() +
   geom_line(aes(y = hgt_filt), color = "blue") +
   geom_line(aes(y = hgt_base*env), color = "red") +
   geom_line(aes(y = env), color = "red")

SH[, mean(env), by = .(lon, k)] %>% 
   ggplot(aes(lon, V1)) + 
   geom_line(aes(color = factor(k)))

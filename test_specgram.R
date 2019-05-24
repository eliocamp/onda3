res <- 2.5
t <-  seq(0, 360 - res, by = res)

# Modelación de la amplitud (onda 1 y onda 2)
amp <- 1.5 + 0.4*cos(t*pi/180 - 30*pi/180)  + 0.6*cos(2*t*pi/180 - 25*pi/180) 


# Onda 3 pura
carrier <- cos(3*(t - 10)*pi/180) + 0.8*cos(5*t*pi/180)

# Señal total 
x <- carrier*amp  

f <- signal::specgram(x, n = length(x)/3, Fs = res, overlap = length(x)/3 - 1)

f$f <- length(x)*(f$f)/res

by =  diff(range(t))/(length(f$t))
ts <- seq(t[1] + by/2, t[length(t)] - by/2, by = by)
f$t <- ts
dimnames(f$S) <- list(f = f$f, t = f$t)



f.df <- as.data.table(melt(abs(f$S)/nrow(f$S)*2))


ggplot(f.df[f == 3], aes(t, value)) +
   geom_line(data = data.frame(t = t, value = x), color = "gray") +
   # geom_line() +
   # geom_point()  +
   geom_line(data = data.frame(t = t, value = amp), color = "red") +
   geom_line(data = data.frame(t = t, value = WaveEnvelope(x))) +
   geom_line(data = data.frame(t = t, value = WaveEnvelope(x, 1:5))) +
   # geom_line(data = data.frame(t = t, value = PeriodicWavelet(x, 3)[[1]]))  +
   NULL


fortify.list <- function(model, data, ... ){
   as.data.frame(model)
}

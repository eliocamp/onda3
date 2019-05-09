OPP <- function(data, Tau) {
  C <- function(tau, g) {
    cov(g[(1 + tau):nrow(g), ], 
        g[1:(nrow(g) - tau), ])
  }
  
  M <- function(Tau = 10, g){
    out <- C(0, g)
    for (t in seq_len(Tau)) {
      out <- out + C(t, g)
    }
    return(out)
  }
  
  C0 <- C(0, data)
  M_ <- M(Tau, data)
  
  X <- solve(C0) %*% (M_ %*% t(M_))
  eig <- eigen(X)
  
  data_prime <- as.data.frame(data %*% eig$vectors)
  colnames(data_prime) <- paste0("V", seq_len(ncol(data_prime)))
  
  r <- C0 %*% eig$vectors
  
  return(list(T1 = eig$values,
              data = data_prime,
              OPP = r))
}


params <- list(sigma = 10, b = 8/3, r = 28)
y <- c(x = 0.1, y = 1, z = 1)


Lorenz <- function(t, y, parms,...) {
  x_dot <- params$sigma*(y[2] - y[1])
  y_dot <- -y[3]*y[1] + params$r*y[1] - y[2]
  z_dot <- y[1]*y[2] - params$b*y[3]
  
  list(c(x_dot, y_dot, z_dot))
}


sol <- deSolve::rk4(y, seq(0, 1800, by = 1/24), Lorenz, params)

g <- as.matrix(sol)[, -1]
g <- scale(g, scale = FALSE)


lorenz_opp <- OPP(g, 100)

t1 <- as.data.frame(lorenz_opp$data)

colnames(t1) <- c("x", "y", "z")
library(plotly)

plotly::plot_ly(t1, x = ~x, y = ~y, z = ~z) %>% 
  plotly::add_paths()

library(metR)
library(data.table)
gh <- ReadNetCDF("~/DATOS/NCEP Reanalysis/hgt.mon.mean.nc", c(gh = "hgt"),
                 subset = list(level = 500, lat = 20:80,
                               time = c("1950-01-01", "1999-12-31")))

gh[, gh.w := Anomaly(gh), by = .(lon, lat, month(time))][, gh.w := gh.w*sqrt(cos(lat*pi/180))]


eof <- EOF(gh.w ~ lon + lat | time, n = 1:10, data = gh)

opp <- OPP(as.matrix(dcast(eof$right, time ~ PC, value.var = "gh.w")[, -1]), 60)

plot(opp$data$V1, type = "l")

acf(opp$data$V1, 100)



comb <- melt(opp$OPP[, 1])
comb$PC <- rownames(comb)

patt <- as.data.table(comb)[eof$left, on = "PC"]


pat <- patt[, .(T1 = sum(gh.w*value)), by = .(lon = ConvertLongitude(lon), lat)]


pat <- ggperiodic::periodic(pat, lon = c(-180, 180))

world <- subset(fortify(rnaturalearth::ne_coastline()), 
                lat > 0)




ggplot(pat, aes(lon, lat))  +
  geom_contour_fill(aes(z = T1)) +
  geom_contour2(aes(z = T1), size = 0.2) +
  geom_path(data = world, aes(long, lat, group = group), color = "gray20") +
  scale_y_latitude(limits = c(20, 90)) +
  scale_fill_divergent() +
  ggalt::coord_proj("+proj=laea +lat_0=90 +lon_0=-180")

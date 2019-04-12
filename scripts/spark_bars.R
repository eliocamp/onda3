spark_bars <- function(x, midpoint = 0, colors = FALSE) {
   lookup_bars <- c("0" = "\033[4m \033[24m",
                    "1" = "\U2582",
                    "2" = "\U2583",
                    "3" = "\U2584",
                    "4" = "\U2585",
                    "5" = "\U2586",
                    "6" = "\U2587")

   normalised <- rep(0, length(x))
   M <- max(x)
   m <- min(x)
   x <- x - midpoint
   
   normalised[x != 0] <- as.numeric(cut(abs(x[x != 0]), breaks = 6))
   normalised[x < 0] <- 7 - normalised[x < 0]
   
   spark <- rep(" ", length = length(normalised))
   
   spark <- lookup_bars[as.character(normalised)]
   
   spark[x < 0] <- paste0("\033[7m", spark[x < 0], "\033[27m")  # inverse
      
      # crayon::inverse(spark[x < 0])
   
   if (isTRUE(colors)) {
      spark[x > 0] <- paste0("\033[31m", spark[x > 0], "\033[39m")  # red
      spark[x < 0] <- paste0("\033[34m", spark[x < 0], "\033[39m")  # blue
   }
   
   
   print_spark(spark, x)
   return(invisible(x))
}


print_spark <- function(spark, x) {
   width <- options()$width
   N <- length(spark)
   start <- seq(1, N, by = width)
   end <- c(start - 1, N)[-1]
   for (s in seq_along(start)) {
      chunk <- seq(start[s], end[s])
      this_spark <- spark[chunk]
      this_x <- x[chunk]
      print_spark_oneline(this_spark, this_x)
      if (length(start) > 1 & s < length(start)) cat("\n")
   }
}

print_spark_oneline <- function(spark, x) {
   # Print positives
   for (i in seq_along(spark)) {
      if (x[i] >= 0) {
         cat(spark[i])
      } else {
         cat(" ")
      }
   }
   cat("\n")
   
   # Print negatives
   for (i in seq_along(spark)) {
      if (x[i] < 0) {
         cat(spark[i])
      } else {
         cat(" ")
      }
   }
}

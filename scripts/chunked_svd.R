chunked_svd <- function(fun_x, nu = 1, nv = nu) {
   sub_Ss <- list()
   data <- 1
   chunk <- 1
   
   while (!is.null(data)) {
      message(paste0("Processing chunk ", chunk, "..."))
      message("   Getting chunk...")
      data <- fun_x(chunk)
      if (!is.null(data)) {
         message("   Processing chunk...")
         sub_Ss[[chunk]] <- svd(data)   
         chunk <- chunk + 1
      }
   }
   message("Computing svd...")
   Y <- Reduce(rbind, 
               lapply(sub_Ss, function(S) {
                  diag(S$d) %*% t(S$v)
               }))
   
   svd_Y <- irlba::irlba(Y, nu = nu, nv = nv)
   
   
   U <- Reduce(rbind, lapply(seq_along(sub_Ss), function(i) {
      S <- sub_Ss[[i]]
      n <- ncol(S$u)
      r <- S$u %*% svd_Y$u[seq_len(n), , drop = FALSE]
      
      svd_Y$u <<- svd_Y$u[-seq_len(n), , drop = FALSE]
      r
      
   }))
   
   # Return SVDs for each chunk
   i <- seq_len(nu)
   sub_Ss <- lapply(sub_Ss, function(S) {
      list(
         d = S$d[i],
         u = S$u[, i, drop = FALSE],
         v = S$v[, i, drop = FALSE]
      )
   })
   
   list(d = svd_Y$d,
        v = svd_Y$v,
        u = U,
        chunks = sub_Ss
   )
}




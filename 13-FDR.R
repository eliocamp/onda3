

data[lat == 50 & lon == 150] %>% 
   lm(var ~ year, data = .) -> fit

# FDR
point[, p.val := pt(abs(estimate)/std.error, df, lower.tail = FALSE)]
p.vals <- a$p.val

FDR <- function(p.vals, alpha = 0.05) {
   M <- length(p.vals)
   id <- order(p.vals)
   p.vals_ordered <- p.vals[id]
   
   a <- alpha * seq_along(p.vals)/M
   signif <- p.vals_ordered <= a
   signif[order(id)]
}


point[term == "year"] %>%
   # .[lon %between% c(150, 250) & lat %between% c(25, 50)] %>%
   .[, p.val := pt(-abs(estimate)/std.error, df)*2] -> a
   .[, signif := FDR(p.val)] %>%  
   ggplot(aes(lon, lat)) +
   geom_contour_fill(aes(z = estimate)) +
   # stat_subset(aes(subset = p.val < 0.05), geom = "point", size = 0.1, alpha = 0.3) +
   stat_subset(aes(subset = signif), geom = "point", size = 0.1, alpha = 0.3) +
   geom_world +
   scale_fill_divergent("Trend in U(300hpa)") +
   coord_quickmap() +
   metR:::theme_field() 

load("mse_results.Rdata")
coarse.pi <- results$coarse.pi
coarse.rho <- results$coarse.rho
rm(results)

plot.width <- 6
plot.height <- 7
legend.inset <- 0.03

line.types <- 1:10
line.colors <- c("black", "red", "blue", "orange", "green", "cyan")
line.width <- 2

par (mar = c(3,3,2,1), mgp = c(2, 0.7, 0), tck = -0.01) 


rmse.plot <- function(panel, col.list, relative = TRUE, 
                      legend.pos = "topright"){
  panel[,-c(1:3)] <- apply(panel[,-c(1:3)], c(1,2), sqrt)
  panel <- as.data.frame(panel) #allows $ column selection
  y <- panel[,col.list] 
  if(relative){
    oracle <- pmin(panel$OLS, panel$TSLS)
    y <- y[,!(col.list %in% c("OLS", "TSLS")), drop = FALSE]
    y <- (y - oracle) / oracle * 100
    y.label <- "RMSE Relative to Oracle (\\%)"
  }else{
    y.label <- "RMSE"
  }
  p <- unique(panel$p)
  r <- unique(panel$r)
  if(length(p) != 1){
    x <- p
    x.label <- "$\\pi$"
    fixed.value <- r
    fixed.label <- "\\rho"
  }else{
    x <- r
    x.label <- "$\\rho$"
    fixed.value <- p
    fixed.label <- "\\pi"
  }
  n <- unique(panel$n)
  plot.title <- paste0("$N=", n, ", \\;", 
                       fixed.label, "=", fixed.value, "$")
  
  matplot(x, y, type = 'l', lty = line.types, 
          col = line.colors, lwd = line.width, 
          xlab = x.label, ylab = y.label, 
          main = plot.title)
  legend(legend.pos, bty = "n", names(y), lty = line.types, 
         col = line.colors, lwd = line.width,
         inset = legend.inset)  
}



plot.grid <- function(results, ...){
  nCols <- length(results)
  nRows <- lapply(results, length)[[1]]
  par(mfcol = c(nRows, nCols))
  for(i in 1:nCols){
    for(j in 1:nRows){
      rmse.plot(results[[i]][[j]], ...)
    }
  }
  par(mfcol = c(1,1))
}

baseline <- c("OLS", "TSLS", "FMSC")
relative.DHW <- c("OLS", "TSLS", "FMSC", "DHW90", "DHW95")
relative.all <- c("OLS", "TSLS", "FMSC", "AVG","DHW90", "DHW95")


tikz('RMSE_coarse_pi_baseline.tex',
     width = plot.width, height = plot.height)
  plot.grid(coarse.pi, baseline, relative = FALSE,
            legend.pos = "topleft")
dev.off()

tikz('RMSE_coarse_pi_relative_DHW.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.pi, relative.DHW)
dev.off()

tikz('RMSE_coarse_pi_relative_all.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.pi, relative.all)
dev.off()

tikz('RMSE_coarse_rho_baseline.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, baseline, relative = FALSE)
dev.off()

tikz('RMSE_coarse_rho_relative_DHW.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, relative.DHW, 
          legend.pos = "topleft")
dev.off()

tikz('RMSE_coarse_rho_relative_all.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, relative.all)
dev.off()


#Clean up
rm(coarse.pi, coarse.rho)
rm(line.colors, line.types, line.width)
rm(rmse.plot, plot.grid)
rm(baseline, relative.DHW, relative.all)
rm(plot.width, plot.height)

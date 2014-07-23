load("./Results/mse_results.Rdata")
coarse.pi <- results$coarse.pi
coarse.rho <- results$coarse.rho
rm(results)

#Set parameters defining plot layout, etc.
nRows <- 3
nCols <- 3
plot.width <- 6
plot.height <- 7
legend.inset <- 0.03
line.types <- 1:10
line.colors <- c("black", "red", "blue", "orange", "green", "cyan")
line.width <- 2
par(mar = c(3,3,2,1), mgp = c(2, 0.7, 0), tck = -0.01) 

#A "panel" is a subset of the simulation results in which
#all but one parameter is fixed. Given a panel and a list
#of columns to plot, this function detects which of the
#parameters p and r is being held constant and sets up the
#plot appropriately. The default relative = TRUE gives RMSE
#relative to the oracle, i.e. the infeasible estimator that
#takes the lower RMSE envelope of OLS and TSLS.
rmse.plot <- function(panel, col.list, relative = FALSE, 
                      legend.pos = "topright"){
  
  results.cols <- setdiff(names(panel), c("n", "p", "r"))
  panel[,results.cols] <- sqrt(panel[,results.cols])
  panel[,setdiff(names(panel), c("n", "p", "r"))]
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



#Convert simulation results into lists of "panels" for
#passing to rmse.plot
coarse.pi.params <- expand.grid(n = unique(coarse.pi$n), 
                                p = unique(coarse.pi$p))
coarse.pi <- mapply(function(x, y) 
                    subset(coarse.pi, n == x & p == y),
                    x = coarse.pi.params$n,
                    y = coarse.pi.params$p, 
                    SIMPLIFY = FALSE)

coarse.rho.params <- expand.grid(n = unique(coarse.rho$n), 
                                 r = unique(coarse.rho$r))
coarse.rho <- mapply(function(x, y) 
                     subset(coarse.rho, n == x & r == y),
                     x = coarse.rho.params$n,
                     y = coarse.rho.params$r, 
                     SIMPLIFY = FALSE)


#This function takes a list of panels and repeatedly calls
#rmse.plot to construct a grid of plots. We can pass it 
#any additional arguments we like and they will "fall through"
#to rmse.plot
plot.grid <- function(results, nRows, nCols, ...){
  par(mfrow = c(nRows, nCols))
  lapply(results, rmse.plot, ...)
  par(mfrow = c(1,1))
}


#Columns to include in the different plots
baseline <- c("OLS", "TSLS", "FMSC")
relative.DHW <- c("OLS", "TSLS", "FMSC", "DHW90", "DHW95")
relative.all <- c("OLS", "TSLS", "FMSC", "AVG","DHW90", "DHW95")


#Create plots and export as tikz files
tikz('./Results/RMSE_coarse_pi_baseline.tex',
     width = plot.width, height = plot.height)
  plot.grid(coarse.pi, nRows, nCols, baseline, 
            legend.pos = "topleft")
dev.off()

tikz('./Results/RMSE_coarse_pi_relative_DHW.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.pi, nRows, nCols, relative.DHW)
dev.off()

tikz('./Results/RMSE_coarse_pi_relative_all.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.pi, nRows, nCols, relative.all)
dev.off()

tikz('./Results/RMSE_coarse_rho_baseline.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, nRows, nCols, baseline)
dev.off()

tikz('./Results/RMSE_coarse_rho_relative_DHW.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, nRows, nCols, relative.DHW, 
          legend.pos = "topleft")
dev.off()

tikz('./Results/RMSE_coarse_rho_relative_all.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, nRows, nCols, relative.all)
dev.off()


#Clean up
rm(coarse.pi, coarse.rho, nRows, nCols)
rm(line.colors, line.types, line.width)
rm(rmse.plot, plot.grid)
rm(baseline, relative.DHW, relative.all)
rm(plot.width, plot.height, legend.inset)
rm(coarse.pi.params, coarse.rho.params)
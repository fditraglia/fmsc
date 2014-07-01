load("./Results/mse_results.Rdata")
coarse.gamma <- results$coarse.gamma
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
rmse.plot <- function(panel, col.list, relative = TRUE, 
                      legend.pos = "topright"){
  
  results.cols <- setdiff(names(panel), c("n", "g", "r"))
  panel[,results.cols] <- sqrt(panel[,results.cols])
  y <- panel[,col.list] 
  
  if(relative){
    oracle <- pmin(panel$Valid, panel$Full)
    y <- y[,!(col.list %in% c("Valid", "Full")), drop = FALSE]
    y <- (y - oracle) / oracle * 100
    y.label <- "RMSE Relative to Oracle (\\%)"
  }else{
    y.label <- "RMSE"
  }
  
  g <- unique(panel$g)
  r <- unique(panel$r)
  
  if(length(g) != 1){
    x <- g
    x.label <- "$\\gamma$"
    fixed.value <- r
    fixed.label <- "\\rho"
  }else{
    x <- r
    x.label <- "$\\rho$"
    fixed.value <- g
    fixed.label <- "\\gamma"
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
coarse.gamma.params <- expand.grid(n = unique(coarse.gamma$n), 
                                g = unique(coarse.gamma$g))
coarse.gamma <- mapply(function(x, y) 
                    subset(coarse.gamma, n == x & g == y),
                    x = coarse.gamma.params$n,
                    y = coarse.gamma.params$g, 
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
baseline <- c("Valid", "Full", "FMSC")
rel.pos <- c("FMSC", "posFMSC")
rel.J <- c("FMSC", "J90", "J95")
rel.MSC <- c("FMSC", "AIC", "BIC", "HQ")
rel.combMSC <- c("FMSC", "combAIC", "combBIC", "combHQ")


#Create plots and export as tikz files

#Coarse grid for Gamma
tikz('./Results/RMSE_coarse_gamma_baseline.tex',
     width = plot.width, height = plot.height)
  plot.grid(coarse.gamma, nRows, nCols, baseline, relative = FALSE,
            legend.pos = "topleft")
dev.off()
 
tikz('./Results/RMSE_coarse_gamma_rel_pos.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.gamma, nRows, nCols, rel.pos)
dev.off()

tikz('./Results/RMSE_coarse_gamma_rel_J.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.gamma, nRows, nCols, rel.J)
dev.off()

tikz('./Results/RMSE_coarse_gamma_rel_MSC.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.gamma, nRows, nCols, rel.MSC)
dev.off()

tikz('./Results/RMSE_coarse_gamma_rel_combMSC.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.gamma, nRows, nCols, rel.combMSC)
dev.off()



#Coarse Grid for Rho
tikz('./Results/RMSE_coarse_rho_baseline.tex',
     width = plot.width, height = plot.height)
  plot.grid(coarse.rho, nRows, nCols, baseline, relative = FALSE,
            legend.pos = "topleft")
dev.off()
 
tikz('./Results/RMSE_coarse_rho_rel_pos.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, nRows, nCols, rel.pos)
dev.off()

tikz('./Results/RMSE_coarse_rho_rel_J.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, nRows, nCols, rel.J)
dev.off()

tikz('./Results/RMSE_coarse_rho_rel_MSC.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, nRows, nCols, rel.MSC)
dev.off()

tikz('./Results/RMSE_coarse_rho_rel_combMSC.tex',
     width = plot.width, height = plot.height)
plot.grid(coarse.rho, nRows, nCols, rel.combMSC)
dev.off()

#Clean up
rm(coarse.gamma, coarse.rho, nRows, nCols)
rm(line.colors, line.types, line.width)
rm(rmse.plot, plot.grid)
rm(baseline, rel.pos, rel.J, rel.MSC, rel.combMSC)
rm(plot.width, plot.height, legend.inset)
rm(coarse.gamma.params, coarse.rho.params)
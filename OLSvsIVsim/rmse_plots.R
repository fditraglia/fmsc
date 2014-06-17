#Load and unpack rmse results
load("rmse_results.Rdata")
coarse.pi <- results$coarse.pi
coarse.rho <- results$coarse.rho
rm(results)

line.types <- 1:10
line.colors <- c("black", "red", "blue", "orange", "green", "cyan")
line.width <- 2

op <- par()
par (mar = c(3,3,2,1), mgp = c(2, 0.7, 0), tck = -0.01) 

#Input a panel (matrix)
#and a list of columns we want to plot
#Function figures out which parameter (p or r) is constant
#and draws the plot accordingly
#Relative RMSE as an option.
rmse.plot <- function(panel, col.list, relative = TRUE){
  panel[,-c(1:3)] <- apply(panel[,-c(1:3)], c(1,2), sqrt)
  panel <- as.data.frame(panel) #allows $ column selection
  y <- panel[,col.list] 
  if(relative){
    oracle <- pmin(panel$b.ols, panel$b.tsls)
    y <- y[,!(col.list %in% c("b.ols", "b.tsls")), drop = FALSE]
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
  plot.title <- paste0("$N=", n, ", \\quad ", 
                       fixed.label, "=", fixed.value, "$")
  
  matplot(x, y, type = 'l', lty = line.types, 
          col = line.colors, lwd = line.width, 
          xlab = x.label, ylab = y.label, 
          main = plot.title)
  legend("topright", names(y), lty = line.types, 
         col = line.colors, lwd = line.width)  
}

test.panel <- coarse.pi[[1]][[1]]
test.cols <- c("b.ols", "b.tsls", "b.fmsc", 
               "b.star","b.DHW90", "b.DHW95")

rmse.plot(test.panel, test.cols)
rmse.plot(test.panel, c("b.ols", "b.tsls", "b.fmsc"), relative = FALSE)

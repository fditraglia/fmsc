#Changed January 16th, 2012 to make the figures easier to see for presentations

#Original Version: September 24th, 2011

#This script makes a figure illustrating the non-normality of the distribution of estimators post-selection


#Home computer directory
working.dir <- "/Users/fditraglia/Dropbox/PhD Dissertation/Job Market Paper/Simulations/Linear IV with Coverage Probs/Selection Uncertainty Figure"

setwd(working.dir)

#Load workspace
load("workspace_image.RData")

#Function to plot the density of a given post-selection estimator alongside that of b.full and b.valid

dist.plot <- function(b.estimator, name){
  
  f.valid <- density(b.valid)
  f.full <- density(b.full)
  f.criterion <- density(b.estimator)
  
  xmin <- min(min(f.valid$x), min(f.full$x), min(f.criterion$x))
  xmax <- max(max(f.valid$x), max(f.full$x), max(f.criterion$x))
  ymin <- min(min(f.valid$y), min(f.full$y), min(f.criterion$y))
  ymax <- max(max(f.valid$y), max(f.full$y), max(f.criterion$y))
  
  plot(f.valid$x, f.valid$y, xlim = c(-0.5, xmax), ylim = c(ymin, ymax), xlab = '', ylab = 'Density', type = 'l', lty = 2, col = 'blue', lwd = 3)
  
  points(f.full$x, f.full$y, type = 'l', lwd = 3, lty = 2, col = 'red')

  points(f.criterion$x, f.criterion$y, type = 'l', lty = 1, lwd = 3)

  legend('topleft', c('Valid', 'Full', name), col = c('blue', 'red', 'black'), lty = c(2,2,1), lwd = c(3,3,3))
  
}#END dist.plot



#Make plots and write to file

pdf('GMM_BIC_present.pdf')
dist.plot(b.BIC, 'GMM-BIC')
dev.off()

pdf('GMM_AIC_present.pdf')
dist.plot(b.AIC, 'GMM-AIC')
dev.off()

pdf('GMM_HQ_present.pdf')
dist.plot(b.HQ, 'GMM-HQ')
dev.off()

pdf('FMSC_present.pdf')
dist.plot(b.FMSC, 'FMSC')
dev.off()

pdf('J_90_present.pdf')
dist.plot(b.J.90, 'J-test (90%)')
dev.off()

pdf('J_95_present.pdf')
dist.plot(b.J.95, 'J-test (95%)')
dev.off()

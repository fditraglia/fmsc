# matplot(r.seq, apply(bar, 2, sqrt), xlab = 'Cor(e,v)', ylab = 'RMSE',
#         type = 'l', lty = 1, lwd = 2)
# legend("topleft", c("OLS", "2SLS", "FMSC", "AVG", "DHW90", "DHW95"), fill = 1:6)
# 
# #Plot expressed in percentage points relative to the oracle
# oracle <- pmin(bar[,1], bar[,2])
# rel.rmse <- (bar - oracle)/oracle * 100
# matplot(r.seq, rel.rmse[,-c(1:2)], xlab = 'Cor(e,v)', ylab = 'RMSE Relative to Oracle (%)', type = 'l', lty = 1:4, lwd = 4, col = c("black", "blue", "red", "orange"))
# legend("topright", c("FMSC", "AVG", "DHW90", "DHW95"), fill = c("black", "blue", "red", "orange"), lty = 1:4, lwd = 3)
load("rmse_results.Rdata")
params <- out[,c("g", "r", "n")]
mse <- out[,setdiff(names(out), c("g", "r", "n"))]
rmse <- sqrt(mse)
cbind(params, round(rmse, 2))
cbind(params, diff  = with(rmse, round(Full - Valid, 2)))
cbind(params, FMSC = with(rmse, round(FMSC, 2)))
cbind(params, diff = with(rmse, round(FMSC - Valid, 2)))
cbind(params, diff = with(rmse, round(FMSC - Full, 2)))
oracle <- with(rmse, pmin(Valid, Full))
cbind(params, round(100 * (rmse - oracle)/oracle))

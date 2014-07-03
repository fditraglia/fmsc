w.rule <- c(0, 1, 0)
w.malfal <- c(0, 0, 1)

g <- function(M, w) t(w) %*% M %*% w

rule.sqbias <- apply(fmsc.ingredients$sqbias, 3, g, w.rule)
rule.avar <- apply(fmsc.ingredients$avar, 3, g, w.rule)
fmsc.rule <- rule.sqbias + rule.avar
fmsc.pos.rule <- pmax(rule.sqbias, 0) + rule.avar

malfal.sqbias <- apply(fmsc.ingredients$sqbias, 3, g, w.malfal)
malfal.avar <- apply(fmsc.ingredients$avar, 3, g, w.malfal)
fmsc.malfal <- malfal.sqbias + malfal.avar
fmsc.pos.malfal <- pmax(malfal.sqbias, 0) + malfal.avar

rule <- cbind(FMSC = fmsc.rule, posFMSC = fmsc.pos.rule, 
              est = drop(t(w.rule) %*% fmsc.ingredients$est))

malfal <- cbind(FMSC = fmsc.malfal, posFMSC = fmsc.pos.malfal, 
              est = drop(t(w.malfal) %*% fmsc.ingredients$est)) 

fmsc.values <- list(rule = rule, malfal = malfal)
  

#Clean Up
rm(w.rule, w.malfal, g)
rm(rule.sqbias, rule.avar, fmsc.rule, fmsc.pos.rule)
rm(malfal.sqbias, malfal.avar, fmsc.malfal, fmsc.pos.malfal)
rm(rule, malfal)

#The FMSC values calculated here are slightly different from those in the original version of the paper for some of the candidates. This is because I originally used submatrices of the full variance matrix for *all* candidates except for the valid specification, whereas I now calculate the variance matrix for *each candidate* separately. Results for the full and valid specifications are completely unchanged.
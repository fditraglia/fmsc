TeXtable <- function(Rtable){
  nRow <- nrow(Rtable)
  nCol <- ncol(Rtable)
  body <- rbind(colnames(Rtable), apply(Rtable, c(1,2), format))
  body <- apply(body, c(1,2), function(x) paste0('$', x, '$'))
  rownames(body) <- NULL
  colnames(body) <- NULL
  body <- apply(body, 1, function(x) paste(x, collapse = ' & '))
  header <- body[1]
  body <- body[-1]
  body[nRow %/% 2] <- paste0("$\\pi\\quad$", body[nRow %/% 2])#rowLegend
  body <- paste(body, collapse = "\\\\ \n")
  colLegend <- paste0("&\\multicolumn{", nCol - 1, "}{c}{$\\rho$}") 
  header <- paste("\\hline\\hline\n", colLegend, "\\\\ \n", header, "\\\\ \n \\hline")
  header <- paste0("\\begin{tabular}{r|", paste(rep("r", nCol - 1), collapse = ""), "}\n", header)
  footer <- "\\\\ \n \\hline \n \\end{tabular}"
  return(paste0(header, body, footer))
}

unlist(lapply(cover.list$FMSC_correct, TeXtable))

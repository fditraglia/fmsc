#Frank DiTraglia
#Last Updated: June 14th, 2014

#This script constructs a table of 2SLS Results for all instrument sets considered in my empirical example. It should not be run on its own: it is called by run_empirical_example.R

#Since I need to do this twice (two panels of results so it fits on one page), I've wrapped everything into a function. 
make.table <- function(range){

  #2SLS Results Panel
  results.raw <- do.call("cbind", tsls.summaries[range])
  results <- results.raw
  rNames <- rownames(results)
  nCols <- ncol(results)
  results <- apply(results, c(1,2), format, nsmall = 2) #format allows us to keep trailing zeros
  results <- apply(results, c(1,2), function(x) paste0('$', x, '$'))
  results <- cbind(rNames, results)
  results <- apply(results, 1, function(x) paste(x, collapse = ' & '))
  names(results) <- NULL
  results <- paste(results, collapse = "\\\\ \n")
  
  #Instrument Blocks Panel
  blocks.raw <- instrument.blocks[range]
  full <- unique(unlist(blocks.raw))
  blocks <- lapply(blocks.raw, function(x) replace(rep(NA, length(full)), match(x, full), x))
  blocks <- t(do.call("rbind", blocks))
  blocks[is.na(blocks)] <- ""
  blocks <- apply(blocks, c(1,2), function(x) paste0('\\multicolumn{2}{c}{', x, "}"))
  blocks <- apply(blocks, 1, function(x) paste(x, collapse = " & "))
  blocks <- paste('&', blocks)
  blocks <- paste(blocks, collapse = "\\\\ \n")

  #Head of Table
  panel.numbers <- paste0('\\multicolumn{2}{c}{', range, '}')
  panel.numbers <- paste(panel.numbers, collapse = ' & ')
  panel.numbers <- paste('&', panel.numbers)
  cNames <- colnames(results.raw)
  cNames <- paste0('\\multicolumn{1}{c}{\\emph{', cNames, '}}')
  cNames <- paste(cNames, collapse = ' & ')
  cNames <- paste('&', cNames)
  table.head <- paste0('\\hline \\hline \n', panel.numbers, '\\\\ \n', cNames, '\\\\ \n \\hline \n') 

  out <- paste0(table.head, ' \n', results, '\\\\ \n', blocks, '\\\\ \n \\hline')
  out <- paste0('\\begin{tabular}{l', paste(rep('r', nCols), collapse = ""), '}\n', out)
  out <- paste0(out, '\n\\end{tabular}')
  return(out)
}

#Assemble full table from two separate panels
panel1 <- make.table(1:6)
panel2 <- make.table(7:12)
full.table <- paste(panel1, '\n \\vspace{2em}', panel2)

#Make "MalfalSq" and "RuleSq" look prettier
#gsub("MalfalSq", full.table, "\\emph{malfal}$^2$")
#gsub("RuleSq", full.table, "\\emph{rule}$^2$")

cat(full.table, file = "table_2SLS_results.tex")
rm(full.table, panel1, panel2, make.table)

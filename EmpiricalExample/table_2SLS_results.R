#Frank DiTraglia
#Last Updated: June 14th, 2014

#This script constructs a table of 2SLS Results for all instrument sets considered in my empirical example. It should not be run on its own: it is called by run_empirical_example.R

range <- 1:2
foo.raw <- do.call("cbind", tsls.summaries[range])
foo <- foo.raw
rNames <- rownames(foo)
nCols <- ncol(foo)
foo <- apply(foo, c(1,2), format, nsmall = 2) #format works like as.character but with more precise formatting commands. We need it here to keep the trailing zeros from being deleted
foo <- apply(foo, c(1,2), function(x) paste0('$', x, '$'))
foo <- cbind(rNames, foo)
foo <- apply(foo, 1, function(x) paste(x, collapse = ' & '))
names(foo) <- NULL
foo <- paste(foo, collapse = "\\\\ \n")
cat(foo)

bar.raw <- instrument.blocks[range]
full <- unique(unlist(bar.raw))
bar <- lapply(bar.raw, function(x) replace(rep(NA, length(full)), match(x, full), x))
bar <- t(do.call("rbind", bar))
bar[is.na(bar)] <- ""
bar <- apply(bar, c(1,2), function(x) paste0('\\multicolumn{2}{c}{', x, "}"))
bar <- apply(bar, 1, function(x) paste(x, collapse = " & "))
bar <- paste(bar, collapse = "\\\\ \n")
cat(bar)

#Head of Table
panel.numbers <- paste0('\\multicolumn{2}{c}{', range, '}')
panel.numbers <- paste(panel.numbers, collapse = ' & ')
panel.numbers <- paste('&', panel.numbers)
cNames <- colnames(foo.raw)
cNames <- paste0('\\multicolumn{1}{c}{\\emph{', cNames, '}}')
cNames <- paste(cNames, collapse = ' & ')
table.head <- paste0('\\hline \\hline \n', panel.numbers, '\\\\ \n', cNames, '\\\\ \n \\hline \n') 
cat(table.head)

out <- paste0(table.head, '\\\\ \n', foo, '\\\\ \n', bar, '\\\\ \n \\hline')
out <- paste0('\\begin{tabular}{l', paste(rep('r', nCols), collapse = ""), '}\n', out)
out <- paste0(out, '\n\\end{tabular}')
cat(out)


#There is a very slight discrepancy (a difference of 0.01 in a few places) between these results and those in the original version of my paper. This comes from the way I've chosen to handle missing observations: I now exclude Vietnam from the dataset completely, since we only observe the baseline instruments for this country and I want to hold everything constant except the instruments in this exercise. In constrast CG use Vietnam when it is available, as did I in the original version of the empirical example. Again, the difference in results is miniscule.
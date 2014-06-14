#Frank DiTraglia
#Last Updated: June 14th, 2014

#This script constructs a table of 2SLS Results for all instrument sets considered in my empirical example. It should not be run on its own: it is called by run_empirical_example.R

foo.raw <- tsls.summaries[1:6]
bar.raw <- instrument.blocks[1:6]

foo <- do.call("cbind", foo.raw)
rNames <- rownames(foo)
cNames <- colnames(foo)
nCols <- ncol(foo) + 1
foo <- apply(foo, c(1,2), as.character)
foo <- apply(foo, c(1,2), function(x) paste0('$', x, '$'))
foo <- cbind(rNames, foo)
foo <- rbind( c("", cNames), foo)
foo <- apply(foo, 1, function(x) paste(x, collapse = ' & '))
names(foo) <- NULL
foo <- paste(foo, collapse = "\\\\ \n")
cat(foo)

full <- unique(unlist(bar.raw))
bar <- lapply(bar.raw, function(x) replace(rep(NA, length(full)), match(x, full), x))
bar <- t(do.call("rbind", bar))
bar[is.na(bar)] <- ""
bar <- apply(bar, c(1,2), function(x) paste0('\\multicolumn{2}{c}{', x, "}"))
bar <- apply(bar, 1, function(x) paste(x, collapse = " & "))
bar <- paste(bar, collapse = "\\\\ \n")
cat(bar)

#There is a very slight discrepancy (a difference of 0.01 in a few places) between these results and those in the original version of my paper. This comes from the way I've chosen to handle missing observations: I now exclude Vietnam from the dataset completely, since we only observe the baseline instruments for this country and I want to hold everything constant except the instruments in this exercise. In constrast CG use Vietnam when it is available, as did I in the original version of the empirical example. Again, the difference in results is miniscule.


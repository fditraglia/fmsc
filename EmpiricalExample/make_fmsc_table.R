rule <- round(fmsc.values$rule, 2)
malfal <- round(fmsc.values$malfal, 2)

#Row labels for table of FMSC results
description <-  lapply(instrument.blocks, 
                       function(x) setdiff(x, "Baseline"))
description <- lapply(description, 
                      function(x) paste(x, collapse = ", "))
description <- unlist(description)
description[1] <- "Valid"
description[8] <- "Full"
description <- paste(paste0("(", 1:8, ")"), description)

#Put together the TeX table
body <- cbind(format(malfal), format(rule))
body <- apply(body, c(1,2), function(x) paste0('$', x, '$'))
colNames <- colnames(body)
body <- cbind(description, body)
body <- apply(body, 1, paste, collapse = " & ")
body <- paste(body, collapse = "\\\\ \n")
head1 <- paste("&", "\\multicolumn{3}{c}{$\\mu = malfal$}&",
               "\\multicolumn{3}{c}{$\\mu = rule$}\\\\ \n")
head2 <- paste("&", paste0(colNames, collapse = " & "), "\\\\ \n")
header <- paste("\\begin{tabular}{lcccccc}\n\\hline\\hline\n",
                head1, head2, "\\hline\n")
footer <- "\\\\ \n \\hline\n \\end{tabular}"
tex.table <- paste(header, body, footer)

#Replace "est" with $\widehat{\mu}$
tex.table <- gsub("est", tex.table, replacement = "$\\widehat{\\mu}$", fixed = TRUE)

#Output TeX table
cat(tex.table, file = "./Results/table_fmsc_values.tex")

#Clean up
rm(body, colNames, description, footer)
rm(head1, head2, header)
rm(tex.table)
rm(malfal, rule)
rm(fmsc.values)

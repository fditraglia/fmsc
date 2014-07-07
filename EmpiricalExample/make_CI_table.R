collapseCI <- function(x) paste0("(", paste0(x, collapse = ", "), ")")

formatCIs <- function(CI.panel){
  apply(format(round(CI.panel, 2)), 1, collapseCI)
} 

malfal.panel <- cbind(FMSC = formatCIs(malfal.FMSC), 
                      posFMSC = formatCIs(malfal.posFMSC))

rule.panel <- cbind(FMSC = formatCIs(rule.FMSC), 
                    posFMSC = formatCIs(rule.posFMSC))

rm(malfal.FMSC, malfal.posFMSC)
rm(rule.FMSC, rule.posFMSC)
rm(collapseCI, formatCIs)

CI.table <- cbind(malfal.panel, rule.panel)

rowNames <- row.names(CI.table)
body <- cbind(rowNames, CI.table)
row.names(body) <- colnames(body) <- NULL

body <- apply(body, 1, paste0, collapse = " & ")
body <- paste0(body, collapse = " \\\\ \n ")

colNames <- paste0(c("", colnames(CI.table)), collapse = " & ")
header <- "\\begin{tabular}{lcccc} \n \\hline \\hline \n & \\multicolumn{2}{c}{\\emph{malfal}} & \\multicolumn{2}{c}{\\emph{rule}}\\\\ \n"
header <- paste0(header, colNames, "\\\\ \n \\hline \n")

footer <- "\\\\ \n \\hline \n\\end{tabular}"

TeX.table <- paste0(header, body, footer, collapse = "\\\\ \n")

# Put a diaresis over the "i" in Naive
TeX.table <- gsub("Naive", x = TeX.table, replacement = "Na\\\"{i}ve", fixed = TRUE)

cat(TeX.table, file = "./Results/CIs_EmpiricalExample.tex")

#Clean up
rm(list = ls())
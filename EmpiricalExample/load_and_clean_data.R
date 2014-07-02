#Frank DiTraglia
#Last Updated: June 13th, 2014

#This script should not be run on its own: it is merely one of the steps from run_empirical_example.R

#This script loads and cleans the dataset from Carstensen & Gundlach (2006), henceforth C&G, and appends two columns not used in the original paper: namely rule^2 and malfal^2. 


#-------------------------------------
#     CODE MISSING VALUES AS NA
#-------------------------------------

#The raw dataset uses -999.999 to indicate missingness
CGdata <- read.csv("Carstensen_Gundlach.csv", 
                   header = TRUE, 
                   na.strings = "-999.999", 
                   stringsAsFactors = FALSE)


#-------------------------------------
#  RENAME COLUMNS TO MATCH C&G PAPER
#-------------------------------------

paper.names <- c("rule", "malfal", "exprop", 
                 "lngdpc", "trade", "latitude", 
                 "coast")

dataset.names <- c("kaufman", "mfalrisk", "exprop2", 
                   "lngdpc95", "frarom", "lat", 
                   "landsea")

key <- data.frame(paper.names, dataset.names, 
                  stringsAsFactors = FALSE)

rm(paper.names, dataset.names)

for(i in 1:nrow(key)){
  col <- which(names(CGdata) == key$dataset.names[i])
  names(CGdata)[col] <- key$paper.names[i]
}

rm(key, col, i)

#-------------------------------------
#       DROP UNNEEDED COLUMNS
#-------------------------------------

keep <- c("lngdpc", "rule", "malfal", "maleco", 
          "lnmort", "frost", "humid", "latitude", 
          "eurfrac", "engfrac", "coast", "trade")

keep.colnumbers <- which(names(CGdata) %in% keep)

CGdata <- CGdata[,keep.colnumbers]

rm(keep, keep.colnumbers)

#-------------------------------------
#    DROP ROWS WITH MISSING VALUES
#-------------------------------------

CGdata <- na.omit(CGdata)

#-------------------------------------
# APPEND rule^2 AND malfal^2 COLUMNS
#-------------------------------------

CGdata <- cbind(CGdata, rule.sq = CGdata$rule^2, 
                        malfal.sq = CGdata$malfal^2)

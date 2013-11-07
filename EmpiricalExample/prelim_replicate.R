#Frank DiTraglia
#September 22nd, 2011

#This script attempts to replicate some of the results from Carstensen and Gundlach (2006) using their dataset.

root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Empirical Examples/Carstensen and Gundlach (Development)"

#Read the data. They appear to use -999.999 as an indicator of missing data
setwd(root.dir)
data <- read.csv("Carstensen_Gundlach.csv", header = TRUE, na.strings = "-999.999")

#Change country from factor data to character data
is.character(data$country) #FALSE
data$country <- as.character(data$country)

attach(data)

#There should be 62 countries for which early settler mortality is available. Of these, 14 are removed because they are likely to provide unreliable statistics, leaving 48. Of the 48, two are removed because they are very small countries (less than 1 million inhabitants in 1990), leaving 46. Of the 46, one is removed because it depends heavily on oil production. Where in the dataset are these countries indicated?
sum(!is.na(lnmort)) #46
sum(!is.na(lnmort.1)) #46
sum(!is.na(logmort1)) #45 - This must be the full set of countries considered
sum(!is.na(logmort2)) #35 - Which are excluded here?

sum(lnmort == lnmort.1, na.rm = TRUE) # 46 - So these are identical
#Double-check visually
cbind(country[(!is.na(lnmort))], country[(!is.na(lnmort.1))]) #Yes

#Which countries?
country[(!is.na(lnmort))]
#"Angola"        
#"Argentina"  
#"Australia"  
#"Bangladesh" 
#"Bolivia"   
#"Brazil"   
#"Burkina Fa" 
#"Cameroon"  
#"Canada"    
#"Chile"      
#"Colombia"  
#"Costa Rica" 
#"Dom. Rep."  
#"Ecuador"    
#"El Salvado" 
#"Guatemala"  
#"Guinea"      
#"Honduras"   
#"Hong Kong"  
#"India"     
#"Indonesia" 
#"Cote d'Ivo" 
#"Jamaica"   
#"Kenya"     
#"Madagascar" 
#"Malaysia" 
#"Mauritania" 
#"Mexico"   
#"Morocco"   
#"New Zealan" 
#"Nigeria"   
#"Pakistan"  
#"Panama"    
#"Paraguay"  
#"Peru"        
#"Senegal"    
#"Singapore" 
#"South Afri" 
#"Sri Lanka" 
#"Tanzania"   
#"Trinidad a" 
#"Tunisia"    
#"Uruguay"    
#"USA"           
#"Venezuela"  
#"Vietnam"    


#Which country is removed to get the list down to 45?
country[(is.na(logmort1))&(!is.na(lnmort))] 
# "Mauritania"
#This is strange, since oil was only discovered in Mauritania in 2001...


#Which countries are excluded to cut it down to 35?
country[(!is.na(logmort1))&(is.na(logmort2))]
#"Angola"     
#"Burkina Fa" 
#"Cameroon"   
#"Guinea"     
#"Hong Kong"  
#"Kenya"     
#"Madagascar" 
#"Nigeria"    
#"Tanzania"   
#"Vietnam" 

detach(data)


#Load libraries for instrumental variables and more general GMM estimation
library(sem)
library(gmm)


#Try to replicate the baseline specifications using the availability of logmort1 to get a group of 45 countries
data.45 <- subset(data, (!is.na(logmort1)))

#Baseline Specification 1: lngdpc95 (log GDP/capita in 1995, called lngdpc in the paper) as the dependent variable; kaufman (measure of institutional quality, called rule in the paper), mfalrisk (measure of malaria prevalence, called malfal in the paper), and a constant as the regressors. The instruments are logmort1 (log early settler mortality), and maleco (a measure of malaria ecology). Constants are included.
baseline1 <- tsls(lngdpc95 ~ kaufman + mfalrisk, ~ logmort1 + maleco, data = data.45)

#Baseline Specification 2: Same as 1 except malrisk replaces malfal.
baseline2 <- tsls(lngdpc95 ~ kaufman + malrisk, ~ logmort1 + maleco, data = data.45)

baseline3 <- tsls(lngdpc95 ~ gadp + mfalrisk, ~ logmort1 + maleco, data = data.45)

baseline4 <- tsls(lngdpc95 ~ exprop2 + mfalrisk, ~ logmort1 + maleco, data = data.45)
  
summary(baseline4)


#Ah wait, these seem to do the trick... But they include one country too many... 
data.46 <- subset(data, (!is.na(lnmort)))

baseline1 <- tsls(lngdpc95 ~ kaufman + mfalrisk, ~ lnmort + maleco, data = data.46)

baseline2 <- tsls(lngdpc95 ~ kaufman + malrisk, ~ lnmort + maleco, data = data.46)

baseline3 <- tsls(lngdpc95 ~ gadp + mfalrisk, ~ lnmort + maleco, data = data.46) 

baseline4 <- tsls(lngdpc95 ~ exprop2 + mfalrisk, ~ lnmort + maleco, data = data.46) 



#According to Kai Carstensen, Mauritania should be the country to remove as it lacks a value for the variable kaufman. See if this is correct...
data$country[is.na(data$kaufman)] #Yup an observation is missing!

#I think the default for 2SLS is to remove missing cases, so when I carried out estimating using data.46 above, this must have automatically removed Mauritania! Test this out by re-running the above models without Mauritania and seeing whether this makes

drop.mauritania <- subset(data.46, country != 'Mauritania')

baseline1 <- tsls(lngdpc95 ~ kaufman + mfalrisk, ~ lnmort + maleco, data = drop.mauritania)

baseline2 <- tsls(lngdpc95 ~ kaufman + malrisk, ~ lnmort + maleco, data = drop.mauritania)

baseline3 <- tsls(lngdpc95 ~ gadp + mfalrisk, ~ lnmort + maleco, data = drop.mauritania) 

baseline4 <- tsls(lngdpc95 ~ exprop2 + mfalrisk, ~ lnmort + maleco, data = drop.mauritania) #This one only gives 41 degrees of freedom, so an observation is missing. Which is it?

#Interesting. It appears that for baseline models 3, which does not use the Kaufman measure, they must be including Mauritania. Is another country missing for the variable gadp?
drop.mauritania$country[is.na(drop.mauritania$gadp)] #Vietnam.

#So I think I have this all figured out. Write a new script that explains everything clearly, and re-labels variables to match the paper.







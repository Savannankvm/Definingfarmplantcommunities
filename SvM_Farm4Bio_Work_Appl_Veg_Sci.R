#S van Mesdag
#05 06 2026
#This R script shows the work on the Farm4Bio dataset, prior to integration with the FSE dataset.

#install packages if need be:

install.packages(vegan)
install.packages(tidyverse)

#load packages

library(vegan)
library(tidyverse)

#Loading the data
#This file has species as columns and each row represents a quadrat sample

Farm4Bio <- read.csv("Farm4Bio_Combined_Standardised_Data_Vo3.csv",  
                       header = TRUE)

print(Farm4Bio)

#Selecting the relevant variables for the CCA

Farm4Bio_variables <- Farm4Bio %>% select(ID:REGION)

Farm4Bio_plants <- Farm4Bio %>% select(Acer.campestre:Zea.mays)

#Running the CCA

Farm4Bio.1.cca <- cca(Farm4Bio_plants ~ FEATURE_CODE, data = Farm4Bio_variables)
Farm4Bio.1.cca

#Main results copied and pasted:

#Call:
#  cca(formula = Farm4Bio_plants ~ FEATURE_CODE, data = Farm4Bio_variables) 

#Partitioning of scaled Chi-square:
#  Inertia Proportion
#Total         12.4780    1.00000
#Constrained    0.1714    0.01373
#Unconstrained 12.3067    0.98627


summary(Farm4Bio.1.cca)
anova.cca(Farm4Bio.1.cca)

#Main results copied and pasted:
#Permutation test for cca under reduced model
#Permutation: free
#Number of permutations: 999

#Model: cca(formula = Farm4Bio_plants ~ FEATURE_CODE, data = Farm4Bio_variables)
#            Df ChiSquare   F Pr(>F)    
#Model      1    0.1714 4.6786  0.001 ***
#  Residual 336   12.3067                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

plot(Farm4Bio.1.cca)

# We are interested in further CCA results as well as the results from the 
#summary and anova. To extract variance explained from the CCA object,
#either the constrained or unconstrained flags in the 'model' term can be used
eig_vals_Farm4Bio <- eigenvals(Farm4Bio.1.cca, model = "constrained")

eig_vals_Farm4Bio

# Now converting them to percentages...
eig_perc_Farm4Bio <- eig_vals_Farm4Bio / sum(eig_vals_Farm4Bio) * 100


## Now pulling out the values for the first two axes for later plotting
cca1_perc_Farm4Bio <- eig_perc_Farm4Bio[1]

print(cca1_perc_Farm4Bio)

##Copied and pasted result:
##CCA1 
##100 

#If you have any trouble with this code, please contact Dr Savanna van Mesdag using
#either her Rothamsted email address, savanna.van-mesdag@rothamsted.ac.uk, or
#her current research profile, such as Researchgate.


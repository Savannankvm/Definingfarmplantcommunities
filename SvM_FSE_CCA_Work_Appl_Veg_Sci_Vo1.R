#S van Mesdag, 04/06/2026

#This code relates to the CCA (vegan) work done to demonstrate the statistically significant
#differences between the Boundary, Verge and Margin habitats on the in-field edge.

#This file is written up for use by anyone who has RStudio installed. Please make sure you have
#a recent version of RStudio installed before running the code in this file to generate results.

##Set Directory.

##Installing the packages if need be:

install.packages("vegan")
install.packages("tidyverse")
install.packages("ggplot2")

library(vegan)
library(tidyverse)
library(ggplot2)

#Loading the plant data. Each row represents a quadrat recorded on a field, in the Boundary,
#Margin or Verge. Taxonomic names were standardised using the species matching function on gbif:
#https://www.gbif.org/tools/species-lookup in June 2024. Only the most common species are included
#in this dataset to refine the Canonical Correspondence Analysis - species that were recorded in 
#>=98%  of the quadrats in the original vegetation dataset were selected for this dataset.
#This meant that plant species that appeared in 82 or fewer quadrats were removed from the
#dataset for the CCA analysis. Quadrats that contained only these removed (<=2%) species were also removed
#from the dataset for the CCA analysis.

#Loading the data:
#This file contains species data as columns, each row is a quadrat sample

BMV_Plants <- read.csv("FSE_BMV_4128_CCA_Plants_INPUT_File_Vo1.csv",  
                       header = TRUE, 
                       colClasses = c("numeric"))

#Loading the associated data, including the data for habitat assignment, Boundary, Verge
#or Margin.
#Each row (other than the header) is a quadrat sample

BMV_Variables <- read.csv("FSE_BMV_4128_CCA_Env_INPUT_File_Vo1.csv",  
                          header = TRUE)

#Running the canonical correspondence analysis, investigating the feature_code, i.e habitat...

cca <- cca(BMV_Plants ~ FEATURE_CODE, data = BMV_Variables)

cca

##Results copied and pasted from console:
##Call: cca(formula = BMV_Plants ~ FEATURE_CODE, data = BMV_Variables)

##Inertia Proportion Rank
##Total         23.6611     1.0000     
##Constrained    0.7760     0.0328    2
##Unconstrained 22.8851     0.9672   79

##Inertia is scaled Chi-square

##Eigenvalues for constrained axes:
##  CCA1   CCA2 
##0.5369 0.2391 

##Eigenvalues for unconstrained axes:
##  CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
##0.6234 0.6113 0.6060 0.5757 0.5728 0.5226 0.4946 0.4753 
##(Showing 8 of 79 unconstrained eigenvalues)


#------
#Now running an anova to look at the statistical significance of the analysis, are the plants
#statistically significantly associated with Boundary, Margin and Verge?

anova.cca(cca)

##Results copied and pasted and neatened from console:
##Permutation test for cca under reduced model
##Permutation: free
##Number of permutations: 999

##Model: cca(formula = BMV_Plants ~ FEATURE_CODE, data = BMV_Variables)
##             Df ChiSquare      F Pr(>F)    
##Model        2       0.776 69.918 0.001 ***
##Residual  4124      22.885                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#For clarification, this anova revealed a high F statistic of 69.918 and a low p value of less than
#0.001, which demonstrates that there is a statistically significant association between
#habitat and plant species presence.

## In order to examine the results of the CCA in more detail, it is useful
#to combine the variables and the plant data together,
#making site IDs explicit for tidy joins later on:

species_data <- BMV_Plants %>%
  as.data.frame() %>%
  rownames_to_column(var = "site_id")

env_data <- BMV_Variables %>%
  as.data.frame() %>%
  rownames_to_column(var = "site_id")

## To check that this has worked, the following code can be used:
glimpse(species_data)
glimpse(env_data)


#It is useful to extract scores from the CCA to examine the results further...
# To obtain site scores (ordination of samples), the following code can be used:
site_scores <- scores(cca, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column(var = "site_id")

# To obtain species scores (to examine the ordination of taxa), the following
#code can be used:
species_scores <- scores(cca, display = "species") %>%
  as.data.frame() %>%
  rownames_to_column(var = "species")

#To obtain environmental (biplot) scores, the following code can be used:
env_scores <- scores(cca, display = "bp") %>%
  as.data.frame() %>%
  rownames_to_column(var = "variable")

# Joining site scores to environmental metadata

site_scores_joined <- site_scores %>%
  left_join(env_data, by = "site_id")

# We are interested in further CCA results as well as the results from the 
#summary and anova. To extract variance explained from the CCA object,
#either the constrained or unconstrained flags in the 'model' term can be used
eig_vals <- eigenvals(cca, model = "constrained")

eig_vals

#Copied and pasted results:
#CCA1      CCA2 
#0.5369134 0.2390657 

# Now converting them to percentages...
eig_perc <- eig_vals / sum(eig_vals) * 100


## Now pulling out the values for the first two axes for later plotting
cca1_perc <- eig_perc[1]
cca2_perc <- eig_perc[2]

print(cca1_perc)

##Copied and pasted result:
##CCA1 
##69.502 

print (cca2_perc)

##Copied and pasted result:
##  CCA2 
##30.498 


#We now want to generate visually pleasing and informative plots to show the
#data. While there are automatic plotting tools in vegan, we will use the plotting
#tools in ggplot2 instead.


# In order to create a more visually pleasing plot,
# We can set up a small variation on the default ggplot theme,
# removing gridlines and adding solid axis lines for all following plots:
theme_neat = theme(panel.grid = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.ticks = element_line(colour = "black"),
                   axis.ticks.length = unit(0.25, "cm")
)

#In order to decide which species to include in the plot, 
#we can look at some of the most 'extreme' species in the CCA plot

options(max.print=1000000)

cca$CCA$v

#In order to include specific species in the plot, we can create a 
#vector for our chosen species...

fig_sp <- c("Crataegus_monogyna", "Veronica_persica",
            "Chaerophyllum_temulum", "Galium_aparine")

#Now to create a plot showing the plant data in relation to Boundary, Margin
#and Verge...

ggplot(data = site_scores_joined, aes(x = CCA1, y = CCA2, colour = FEATURE_CODE)) +
  
  ## Site scores coloured by management type
  geom_point(size = 1, alpha = 0.3) +
  
  ## Habitat centroids
  geom_point(
    data = habitat_centroids,
    aes(
      x = CCA1,
      y = CCA2,
    ),
    shape = 4, size = 5,
    fill = "white", colour = "black", stroke = 1.2) +
  
  ## Species points
  geom_point(
    data = species_scores,
    aes(
      x = CCA1,
      y = CCA2
    ),
    colour = "grey40", size=1, shape=17
  )+
  ## species subset for labelling
  geom_point(
    data = species_scores[which(species_scores$species %in% fig_sp), ],
    aes(
      x = CCA1,
      y = CCA2
    ),
    colour = "grey40", size=3, shape=17
  )+
  
  ##Placing the labels in specific parts of the plot
  annotate("text",
           x = species_scores$CCA1[which(species_scores$species %in% fig_sp)],
           y = species_scores$CCA2[which(species_scores$species %in% fig_sp)],
           label = c("V. persica","C. monogyna","G. aparine", "C. temulum" ),
           fontface="italic",
           hjust=1.1,
           size = 10/.pt,
           colour="grey40")+
  
  
  #Labels for Boundary, Verge and Margin
  annotate("text",
           x = habitat_centroids$CCA1,
           y = habitat_centroids$CCA2,
           label = c("Boundary","Margin","Verge"),
           fontface="bold",
           hjust=-0.25,
           size = 12/.pt)+
  
  ## Plot styling
  theme_minimal() +
  theme_neat +
  
  ## Configure the headings and axis labels
  labs(x = paste0("CCA Axis 1 (", round(cca1_perc, 1), "%)"),
    y = paste0("CCA Axis 2 (", round(cca2_perc, 1), "%)"),
    colour = "Habitat"
  )


#In addition to this plot, we want to look at the vegclust communities.
#Beforheand, the following code can be used to remove the N values.

glimpse(site_scores_joined)

which(site_scores_joined$VEGCLUST_m1.1_dn0.90 == "N")
site_scores_joined$VEGCLUST_m1.1_dn0.90[which(site_scores_joined$VEGCLUST_m1.1_dn0.90 == "N")]  <- 0

site_scores_joined_NAs_removed <- 
  site_scores_joined[site_scores_joined$VEGCLUST_m1.1_dn0.90 != 0, ]

#Now to produce a plot showing the CCA and the vegclust communities...

ggplot(data = site_scores_joined_NAs_removed, aes(x = CCA1, y = CCA2, colour = VEGCLUST_m1.1_dn0.90)) +
  ## Site scores coloured by cluster
  geom_point(size = 1, alpha = 0.3) +
  ## Habitat centroids
  geom_point(
    data = habitat_centroids,
    aes(
      x = CCA1,
      y = CCA2,
    ),
    shape = 4, size = 5,
    fill = "white", colour = "black", stroke = 1.2) +
  ## Species points
  geom_point(
    data = species_scores,
    aes(
      x = CCA1,
      y = CCA2
    ),
    colour = "grey40", size=1, shape=17
  )+
  ## species subset for labelling
  geom_point(
    data = species_scores[which(species_scores$species %in% fig_sp), ],
    aes(
      x = CCA1,
      y = CCA2
    ),
    colour = "grey40", size=3, shape=17
  )+
  
  ##Placing the labels in specific parts of the plot
  annotate("text",
           x = species_scores$CCA1[which(species_scores$species %in% fig_sp)],
           y = species_scores$CCA2[which(species_scores$species %in% fig_sp)],
           label = c("V. persica","C. monogyna","G. aparine", "C. temulum" ),
           fontface="italic",
           hjust=1.1,
           size = 10/.pt,
           colour="grey40")+
  annotate("text",
           x = habitat_centroids$CCA1,
           y = habitat_centroids$CCA2,
           label = c("Boundary","Margin","Verge"),
           fontface="bold",
           hjust=-0.25,
           size = 12/.pt)+
  ## Plot styling
  theme_minimal() +
  theme_neat +
  ## Configure the headings and axis labels
  labs(
    #  title = "Canonical Correspondence Analysis (CCA)",
    #  subtitle = "Arable off-crop cluster communities in Boundary, Margin and Verge",
    x = paste0("CCA Axis 1 (", round(cca1_perc, 1), "%)"),
    y = paste0("CCA Axis 2 (", round(cca2_perc, 1), "%)"),
    colour = "Community"
  )



#We have now produced two plots, one clearly showing the Boundary, Verge and 
#Margin, and the other showing the communities determined through analysis
#using the Vegclust package, which demonstrate some of our key findings.

#If you have any issues running this code or have any further questions regarding the code run
#here, please contact Dr Savanna van Mesdag using the email address 
#savanna.van-mesdag@rothamsted.ac.uk, or contact her through her current
#research profile, such as Researchgate.


#04/06/2026
#Written by S van Mesdag and D Comont
################################################################################
## FSE Vegclust Analysis - Revised Script
## Noise clustering of arable field boundary/margin/verge plant community data
## Indval analysis of plant community data
## Biodiversity index calculation for plant community data


#VEGCLUST ANALYSIS, INCLUDING PARAMATERISATION
################################################################################
## Key points:
##   - dnoise explored within a valid range for chord-normalised data [0, sqrt(2)]
##   - m used as an additional parameter for increasing quadrat assignment rate
##   - Systematic grid search replaces ad hoc parameter changes
##   - IndVal and frequency class tables outputted
##   - FSE clusters applied to the Farm4Bio dataset
############################################################################

#Install packages if need be:

install.packages("indicspecies")
install.packages("patchwork")
install.packages("viridis")
install.packages("vegclust")
install.packages("vegan")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("parallel")
install.packages("patchwork")
install.packages("hillR")

## Load some packages
library(vegclust)
library(vegan)
library(indicspecies)
library(tidyverse)
library(dplyr)
library(parallel)
library(patchwork)
library(viridis)
library(hillR)

################################################################################
## 1. DATA LOADING AND PREPARATION
################################################################################


#############################
## 1. LOAD IN THE DATA

#Set the working directory...

## Read in the data for analysis...
#(This is a wide dataset with species as columns and sites as rows)
FSEData <- read.csv("FSEDataWideUnique_Vo1.csv", header = TRUE)

# Retain species columns only
FSEData <- FSEData %>% select(Chenopodium_album:Ranunculus_ficaria)

# Chord normalisation: constrains all Euclidean distances to [0, sqrt(2) ~ 1.414]
# dnoise must be interpreted relative to this range throughout
FSEchord <- decostand(FSEData, "normalize")


################################################################################
## 2. GRID SEARCH: SYSTEMATIC PARAMETER TUNING
################################################################################

###PLEASE NOTE -
#Much of the code run here takes a long time and uses much processing power.
#If you wish to examine the final key result in more detail, skip to line 197

##
## This is what the web suggests about 'dnoise' and 'm' parameters:
##
## m    = fuzziness exponent. Lower values (~1.05-1.15) sharpen cluster membership,
##        increasing assignment rate. Higher values (~1.2+) spread membership more
##        diffusely, pushing more quadrats into the noise class.
##
## dnoise = noise distance threshold in the chord-normalised Euclidean space.
##          Valid range is approximately 0.75-0.95 for this data type.
##          Values approaching or exceeding 1.0 consume too much of the distance
##          space, leaving insufficient seed objects for new clusters to form.
##          Do NOT use values from un-normalised studies (e.g. dnoise = 50).
##
## min.size = minimum quadrats per cluster. 5 is a sensible default; 4 is
##            acceptable if small genuine communities are suspected.
################################################################################


## Set the dnoise and m gradients to include in the search
m_values      <- c(1.05, 1.10, 1.15, 1.20)
dnoise_values <- c(0.75, 0.80, 0.85, 0.90)

grid_params <- expand.grid(m = m_values, dnoise = dnoise_values)
n_total     <- nrow(FSEchord)

## Detect available cores, leaving one free for the OS
n_cores <- max(1, detectCores() - 1)
cat("Running grid search across", nrow(grid_params), "parameter combinations on",
    n_cores, "cores\n")

## Set up the parallel computing environment
cl <- makeCluster(n_cores)
clusterExport(cl, varlist = c("FSEchord", "grid_params", "n_total"))   # send data to workers
clusterEvalQ(cl, library(vegclust))               # load package on each worker

## And run the grid search
## Each list element returns the full vegclust object plus summary stats
grid_list <- parLapply(cl, seq_len(nrow(grid_params)), function(i) {
  
  m_i      <- grid_params$m[i]
  dnoise_i <- grid_params$dnoise[i]
  
  fit <- tryCatch(
    incr.vegclust(FSEchord, method = "NC", m = m_i, dnoise = dnoise_i, min.size = 5),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(m = m_i, dnoise = dnoise_i,
                n_clusters = NA, n_assigned = NA, pct_assigned = NA,
                vegclust_obj = NULL))
  }
  
  groups_i     <- defuzzify(fit)$cluster
  n_assigned_i <- sum(groups_i != "N")
  
  list(
    m             = m_i,
    dnoise        = dnoise_i,
    n_clusters    = length(unique(groups_i[groups_i != "N"])),
    n_assigned    = n_assigned_i,
    pct_assigned  = round((n_assigned_i / n_total) * 100, 1),
    vegclust_obj  = fit
  )
})

## Stop the clusters
stopCluster(cl)

## Name each list element for easy retrieval: e.g. "m1.15_dn0.80"
names(grid_list) <- sprintf("m%.2f_dn%.2f", grid_params$m, grid_params$dnoise)

## Summary table (no objects — just the stats)
grid_results <- do.call(rbind, lapply(grid_list, function(x) {
  data.frame(m = x$m, dnoise = x$dnoise,
             n_clusters   = x$n_clusters,
             n_assigned   = x$n_assigned,
             pct_assigned = x$pct_assigned)
}))
rownames(grid_results) <- NULL

## Have a look at the grid search results
print(grid_results)


## Make some figures of the fitting results
p1 = ggplot(data = grid_results, aes(x = dnoise, y = pct_assigned, colour = as.character(m))) +
  geom_point() +
  geom_line() +
  labs(x = "dnoise", y= "Percentage of quadrats assigned", colour = "m") +
  theme_bw()

p2 = ggplot(data = grid_results, aes(x = dnoise, y = n_clusters, colour = as.character(m))) +
  geom_point() +
  geom_line() +
  labs(x = "dnoise", y= "Number of clusters", colour = "m") +
  theme_bw()

p1 / p2


###########################
## Save everything neatly
###########################

## And save as a .csv
write.csv(grid_results, file = "FSE_Vegclust_GridSearch_Results.csv", row.names = FALSE)

## Save the list of fitted models from the grid search
saveRDS(grid_list, file = "grid_list.rds")

## Save the plot
pdf("Cluster_results.pdf", width = 7, height = 10)
p1 / p2
dev.off()

## Save the image!
## This writes out the workspace for later retrieval
#Update directory accordingly
save.image("FSE_Vegclust_Revised_workspace_clean.RData")


###############################################################################
## 3. EXTRACT CHOSEN MODEL
###############################################################################

## Assuming we've closed the session
## And definitely don't want to run through the whole grid-search again!

##########################
## Load the data back in
##########################

#Read the object back and assign it to a name of your choice
Search_results <- readRDS("grid_list.rds")


## Edit m_chosen and dnoise_chosen based on grid search output above.
m_chosen      <- 1.10
dnoise_chosen <- 0.90

## Extract the model
FSE.nc.final = Search_results[["m1.10_dn0.90"]][["vegclust_obj"]]


## Check cluster assignments (N = noise/unassigned)
groups_final <- defuzzify(FSE.nc.final)$cluster
table_final  <- table(groups_final)
print(table_final)

#Results copied and pasted:
#groups_final
#F10 F11 F12 F13 F14 F15 F16 F17 F18 F19  F2 F20 F21 F22 F23 F24 F25  F3  F4  F5 
#48 149 110  72 100 134  81 116 105 137  40 101 137 203 241 415 516  54  67  73 
#F6  F7  F8  F9  M1   N 
#44  65  85 100  23 936 

n_assigned_final   <- sum(groups_final != "N")
pct_assigned_final <- round((n_assigned_final / n_total) * 100, 1)
cat("Quadrats assigned:", n_assigned_final, "/", n_total,
    "(", pct_assigned_final, "%)\n")

### Result:
#Quadrats assigned: 3216 / 4152 ( 77.5 %)

## Save as a csv
write.csv(groups_final, file = "FSE_Vegclust_FinalGroups.csv")

## Fuzzy membership matrix (rows = quadrats, cols = clusters + noise)
memb_matrix <- FSE.nc.final$memb

## Internal cluster quality
print(clustvar(FSE.nc.final))                          # fuzzy
print(clustvar(FSE.nc.final, defuzzify = TRUE))        # hard

#Results copied and pasted for hard communities, which best represent the data:
#hard
#M1        F2        F3        F4        F5        F6        F7        F8 
#0.2097465 0.4132858 0.5412683 0.3764911 0.5843558 0.3300898 0.4918133 0.5159265 
#F9       F10       F11       F12       F13       F14       F15       F16 
#0.5856420 0.3034367 0.6240383 0.6017263 0.3888887 0.4344904 0.6026606 0.3756002 
#F17       F18       F19       F20       F21       F22       F23       F24 
#0.4065721 0.4430642 0.4158551 0.7493439 0.2839708 0.4575164 0.3207598 0.1782127 
#F25 
#0.2521938 


################################################################################
## 4. SPECIES CHARACTERISATION
################################################################################
## Using approaches to quantify the defining species associated 
## with each cluster. Two approach used:
## Indicator species analysis
################################################################################

#Install packages if need be:

install.packages("labdsv")
install.packages("tidyverse")
install.packages("dplyr")

##Then, load the packages...

library(labdsv)
library(dplyr)
library(tidyverse)

#Read in the data
#(This is a wide dataset with species as columns and sites as rows)

BMV_Vegclust <- read.csv("BMV_Combined_Metadata_Vo11.csv",  
                         header = TRUE)

#Examining the file...

head(BMV_Vegclust)

#For the Indval analyses it is important to set up the data accordingly, this will be done in the next 
#few lines of code...

Plant_Variables_with_NA <- BMV_Vegclust %>% select(VEGCLUST_m1.1_dn0.90, Chenopodium_album:Ranunculus_ficaria)
Vegclust_Variables_Plants <- Plant_Variables_with_NA %>% drop_na(VEGCLUST_m1.1_dn0.90)
Vegclust_Plants <- Plant_Variables_with_NA %>%select(Chenopodium_album:Ranunculus_ficaria)
Plant_Variables <- Plant_Variables_with_NA %>% select(VEGCLUST_m1.1_dn0.90)

#Making sure everything looks as it should be...

head(Plant_Variables)
tail(Plant_Variables)

head(Vegclust_Plants)

#Making sure the Vegclust community data are the right mode and class...

sapply(Plant_Variables, mode)
sapply(Plant_Variables, class)

Plant_Variables <- data.frame(Plant_Variables)

Vegclust_Plants <- data.frame(Vegclust_Plants)

#These data transformations convert the community categories (i.e F2, F3 etc) to numbers, so that 
#the Indval analysis can process them correctly. Please note that "N", numbered "12", represents
#plants that were unassigned, so they will be ignored in the final results.

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F2")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F2")] <- 2

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F3")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F3")] <- 3

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F4")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F4")] <- 4

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F5")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F5")] <- 5

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F6")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F6")] <- 6

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F7")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F7")] <- 7

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F8")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F8")] <- 8

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F9")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F9")] <- 9

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F10")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F10")] <- 10

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F11")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F11")] <- 11

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F12")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F12")] <- 12

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F13")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F13")] <- 13

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F14")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F14")] <- 14

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F15")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F15")] <- 15

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F16")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F16")] <- 16

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F17")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F17")] <- 17

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F18")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F18")] <- 18

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F19")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F19")] <- 19

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F20")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F20")] <- 20

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F21")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F21")] <- 21

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F22")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F22")] <- 22

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F23")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F23")] <- 23

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F24")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F24")] <- 24

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F25")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "F25")] <- 25

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "M1")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "M1")] <- 26

which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "N")
Plant_Variables$VEGCLUST_m1.1_dn0.90[which(Plant_Variables$VEGCLUST_m1.1_dn0.90 == "N")] <- 27

Plant_Variables$VEGCLUST_m1.1_dn0.90 <- as.numeric(Plant_Variables$VEGCLUST_m1.1_dn0.90)

VEGCLUST <- Plant_Variables %>% select(VEGCLUST_m1.1_dn0.90)


#Making sure that any issues with plant species rows can be resolved, the following lines of code
#are used...
#Running this code using tidyverse will create a new version of the Vegclust_Plants, 
#one that excludes all of the columns that have a sum value of less than 1.

library(tidyverse)
df <- tibble(a = 1:5, b = -1:3, c = 0)

print(df)

selection_criteria <- (colSums(abs(df)) != 0)
df[selection_criteria]

Plant_Variables <- tibble("Vegclust_Plants")

selection_criteria <- (colSums(abs(Vegclust_Plants)) != 0)

print(selection_criteria)

Vegclust_Plants[selection_criteria]

#Now, to run the Indval analyses to look at significant results only...

clustVEGCLUST <- cut(VEGCLUST$VEGCLUST_m1.1_dn0.90, 27, labels=FALSE)
VEGCLUSTIndval <- indval(as.matrix(Vegclust_Plants[selection_criteria]), clustVEGCLUST + clustVEGCLUST, type= "long")
summary(VEGCLUSTIndval)
gr <- VEGCLUSTIndval$maxcls[VEGCLUSTIndval$pval<=0.05]
iv <- VEGCLUSTIndval$indcls[VEGCLUSTIndval$pval<=0.05]
pv <- VEGCLUSTIndval$pval[VEGCLUSTIndval$pval<=0.05]
fr <- apply(Vegclust_Plants[selection_criteria][,-1]>0, 2, sum)[VEGCLUSTIndval$pval<=0.05]
VEGCLUSTIndvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
VEGCLUSTIndvalsummary <- VEGCLUSTIndvalsummary[order(VEGCLUSTIndvalsummary$group, 
                                                     -VEGCLUSTIndvalsummary$indval),]
print(VEGCLUSTIndvalsummary)

##If the code has not been re-run, the results have been copied and pasted below within this 
#R file for easier access...
#In summary, the Indval analyses show Indval or fidelity scores for one or more plant species
#in each community between 1 and 0. The closer the Indval score is to 1, the more that plant
#species is associated with that community. P values in the second representation of the results
#may need to be interpreted with caution, as the communities were already assigned using a Vegclust
#analysis prior to the data being analysed using this labdsv package. It makes sense that the presence
#of certain plant species was not 'down to chance'!

# cluster indicator_value probability
#Sisymbrium_officinale                           1          0.8698       0.001
#Tripleurospermum_maritimum                      1          0.0958       0.001
#Geranium_dissectum                              1          0.0592       0.007
#Lapsana_communis                                1          0.0575       0.007
#Geum_urbanum                                    1          0.0567       0.003
#Malva_sylvestris                                1          0.0518       0.002
#Papaver_rhoeas                                  1          0.0400       0.016
#Ballota_nigra                                   1          0.0335       0.003
#Medicago_sativa_subsp_sativa                    1          0.0249       0.017
#Leucanthemum_vulgare                            1          0.0237       0.019
#Stachys_offici0lis                              1          0.0232       0.019
#Descurainia_sophia                              1          0.0216       0.037
#Malva_neglecta                                  1          0.0210       0.027
#Fallopia_convolvulus                            2          0.7224       0.001
#Lysimachia_arvensis_subsp_arvensis              2          0.0835       0.001
#Atriplex_patula                                 2          0.0294       0.050
#Beta_vulgaris                                   2          0.0262       0.029
#Poa_annua                                       3          0.3678       0.001
#Sonchus_asper                                   4          0.5821       0.001
#Alopecurus_myosuroides                          4          0.0727       0.004
#Chaenorhinum_minus                              4          0.0274       0.014
#Sambucus_nigra                                  5          0.8761       0.001
#Hedera_helix                                    5          0.0768       0.001
#Salix_fragilis                                  5          0.0575       0.002
#Corylus_avellana                                5          0.0493       0.002
#Ilex_aquifolium                                 5          0.0290       0.009
#Bryonia_dioica                                  5          0.0233       0.044
#Malva_sp.                                       5          0.0216       0.024
#Solanum_dulcamara                               5          0.0211       0.041
#Viola_arvensis                                  6          0.7523       0.001
#Myosotis_arvensis                               6          0.0607       0.001
#Fumaria_officinalis                             6          0.0445       0.008
#Torilis_japonica                                6          0.0273       0.005
#Linum_usitatissimum                             6          0.0213       0.029
#Ranunculus_repens                               7          0.6660       0.001
#Cerastium_fontanum                              7          0.0677       0.002
#Trifolium_repens                                7          0.0419       0.006
#Deschampsia_cespitosa                           7          0.0266       0.015
#Phleum_pratense                                 7          0.0254       0.006
#Plantago_lanceolata                             7          0.0242       0.038
#Glyceria_declinata                              7          0.0235       0.014
#Veronica_serpyllifolia                          7          0.0216       0.032
#Chenopodium_album                               8          0.6957       0.001
#Urtica_urens                                    8          0.1825       0.001
#Capsella_bursa.pastoris                         8          0.0690       0.004
#Solanum_nigrum                                  8          0.0610       0.001
#Persicaria_maculosa                             8          0.0555       0.004
#Polygonum_aviculare_agg.                        8          0.0502       0.015
#Amaranthus_retroflexus                          8          0.0300       0.010
#Epilobium_hirsutum                              9          0.8902       0.001
#Filipendula_ulmaria                             9          0.2079       0.001
#Vicia_sativa                                    9          0.0872       0.001
#Centaurea_nigra                                 9          0.0714       0.001
#Scabiosa_sp.                                    9          0.0417       0.001
#Stachys_sylvatica                               9          0.0324       0.037
#Scrophularia_auriculata                         9          0.0208       0.044
#Vicia_sativa_subsp._nigra                       9          0.0208       0.038
#Berula_erecta                                   9          0.0208       0.038
#Rorippa_sylvestris                              9          0.0208       0.040
#Leontodon_hispidus                              9          0.0198       0.042
#Petasites_hybridus                              9          0.0198       0.044
#Arrhenatherum_elatius                          10          0.4436       0.001
#Poa_trivialis                                  10          0.2806       0.001
#Dactylis_glomerata                             10          0.1711       0.001
#Bromus_sterilis                                10          0.1640       0.001
#Holcus_lanatus                                 10          0.0523       0.033
#Reseda_luteola                                 10          0.0268       0.017
#Sherardia_arvensis                             10          0.0253       0.018
#Stellaria_media                                11          0.5991       0.001
#Senecio_vulgaris                               11          0.0609       0.005
#Veronica_hederifolia                           11          0.0316       0.004
#Barbarea_vulgaris                              11          0.0183       0.045
#Rumex_obtusifolius                             12          0.5830       0.001
#Rubus_fruticosus_agg.                          13          0.4536       0.001
#Alnus_glutinosa                                13          0.0360       0.002
#Juncus_effusus                                 13          0.0318       0.006
#Angelica_sylvestris                            13          0.0200       0.041
#Veronica_persica                               14          0.4870       0.001
#Lamium_purpureum                               14          0.1278       0.001
#Sinapis_arvensis                               15          0.8979       0.001
#Raphanus_raphanistrum_subsp._raphanistrum      15          0.0503       0.001
#Heracleum_sphondylium                          16          0.4918       0.001
#Conopodium_majus                               17          0.5839       0.001
#Anthriscus_sylvestris                          18          0.6125       0.001
#Pteridium_aquilinum                            20          0.8515       0.001
#Rosa_canina_agg.                               20          0.0815       0.001
#Teucrium_scorodonia                            20          0.0197       0.042
#Cirsium_arvense                                21          0.4934       0.001
#Galium_aparine                                 22          0.3210       0.001
#Crataegus_monogyna                             23          0.6506       0.001
#Dioscorea_communis                             23          0.0273       0.013
#Urtica_dioica                                  24          0.3700       0.001
#Equisetum_arvense                              25          0.9601       0.001
#Prunella_vulgaris                              25          0.0557       0.001
#Senecio_jacobaea                               25          0.0258       0.037
#Epilobium_parviflorum                          25          0.0258       0.021
#Festuca_ovina                                  25          0.0249       0.027
#Cynosurus_cristatus                            25          0.0141       0.047

#Sum of probabilities                 =  210.584

#Sum of Indicator Values              =  21.67

#Sum of Significant Indicator Values  =  19.02

#Number of Significant Indicators     =  98

#Significant Indicator Distribution

#1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 20 21 22 23 24 25 
#13  4  1  3  8  5  8  7 12  7  4  1  4  2  2  1  1  1  3  1  1  2  1  6 

#                                         group     indval pvalue freq
#Sisymbrium_officinale                         1 0.86983704  0.001  509
#Tripleurospermum_maritimum                    1 0.09576667  0.001   84
#Geranium_dissectum                            1 0.05921914  0.007 1502
#Lapsana_communis                              1 0.05748380  0.007   23
#Geum_urbanum                                  1 0.05665991  0.003  149
#Malva_sylvestris                              1 0.05179117  0.002  308
#Papaver_rhoeas                                1 0.03999389  0.016  729
#Ballota_nigra                                 1 0.03346441  0.003  123
#Medicago_sativa_subsp_sativa                  1 0.02485836  0.017    2
#Leucanthemum_vulgare                          1 0.02374607  0.019    1
#Stachys_offici0lis                            1 0.02317352  0.019    3
#Descurainia_sophia                            1 0.02155172  0.037    1
#Malva_neglecta                                1 0.02103960  0.027    1
#Fallopia_convolvulus                          2 0.72235535  0.001  421
#Lysimachia_arvensis_subsp_arvensis            2 0.08354649  0.001  333
#Atriplex_patula                               2 0.02942085  0.050  463
#Beta_vulgaris                                 2 0.02620137  0.029   16
#Poa_annua                                     3 0.36782305  0.001   65
#Sonchus_asper                                 4 0.58212262  0.001  602
#Alopecurus_myosuroides                        4 0.07270357  0.004   17
#Chaenorhinum_minus                            4 0.02739726  0.014    1
#Sambucus_nigra                                5 0.87605051  0.001 1704
#Hedera_helix                                  5 0.07677160  0.001  886
#Salix_fragilis                                5 0.05746934  0.002   16
#Corylus_avellana                              5 0.04926677  0.002   14
#Ilex_aquifolium                               5 0.02897083  0.009   18
#Bryonia_dioica                                5 0.02334910  0.044   84
#Malva_sp.                                     5 0.02158516  0.024    3
#Solanum_dulcamara                             5 0.02107043  0.041    2
#Viola_arvensis                                6 0.75227296  0.001   56
#Myosotis_arvensis                             6 0.06073056  0.001  191
#Fumaria_officinalis                           6 0.04453253  0.008   59
#Torilis_japonica                              6 0.02725215  0.005   51
#Linum_usitatissimum                           6 0.02125510  0.029   17
#Ranunculus_repens                             7 0.66598633  0.001   33
#Cerastium_fontanum                            7 0.06767663  0.002   17
#Trifolium_repens                              7 0.04191245  0.006  215
#Deschampsia_cespitosa                         7 0.02663326  0.015   10
#Phleum_pratense                               7 0.02541343  0.006  392
#Plantago_lanceolata                           7 0.02420319  0.038   71
#Glyceria_declinata                            7 0.02352941  0.014    1
#Veronica_serpyllifolia                        7 0.02161515  0.032    1
#Chenopodium_album                             8 0.69570376  0.001   40
#Urtica_urens                                  8 0.18253618  0.001 1184
#Capsella_bursa.pastoris                       8 0.06900011  0.004   15
#Solanum_nigrum                                8 0.06095434  0.001  129
#Persicaria_maculosa                           8 0.05553618  0.004   94
#Polygonum_aviculare_agg.                      8 0.05020932  0.015  295
#Amaranthus_retroflexus                        8 0.03000000  0.010  175
#Epilobium_hirsutum                            9 0.89022800  0.001  347
#Filipendula_ulmaria                           9 0.20794055  0.001   24
#Vicia_sativa                                  9 0.08718037  0.001  112
#Centaurea_nigra                               9 0.07139602  0.001   16
#Scabiosa_sp.                                  9 0.04166667  0.001    1
#Stachys_sylvatica                             9 0.03239071  0.037   77
#Scrophularia_auriculata                       9 0.02083333  0.044    2
#Vicia_sativa_subsp._nigra                     9 0.02083333  0.038    2
#Berula_erecta                                 9 0.02083333  0.038    1
#Rorippa_sylvestris                            9 0.02083333  0.040   11
#Petasites_hybridus                            9 0.01982907  0.044    1
#Leontodon_hispidus                            9 0.01981707  0.042    2
#Arrhenatherum_elatius                        10 0.44356881  0.001   38
#Poa_trivialis                                10 0.28058787  0.001  540
#Dactylis_glomerata                           10 0.17112861  0.001  475
#Bromus_sterilis                              10 0.16397746  0.001    5
#Holcus_lanatus                               10 0.05228867  0.033  191
#Reseda_luteola                               10 0.02684564  0.017   10
#Sherardia_arvensis                           10 0.02532754  0.018    8
#Stellaria_media                              11 0.59910317  0.001  292
#Senecio_vulgaris                             11 0.06093240  0.005   79
#Veronica_hederifolia                         11 0.03163742  0.004    9
#Barbarea_vulgaris                            11 0.01834710  0.045    1
#Rumex_obtusifolius                           12 0.58301548  0.001  193
#Rubus_fruticosus_agg.                        13 0.45355121  0.001  213
#Alnus_glutinosa                              13 0.03599911  0.002    1
#Juncus_effusus                               13 0.03184366  0.006    3
#Angelica_sylvestris                          13 0.02000000  0.041    5
#Veronica_persica                             14 0.48702047  0.001    9
#Lamium_purpureum                             14 0.12778090  0.001    1
#Sinapis_arvensis                             15 0.89785365  0.001  292
#Raphanus_raphanistrum_subsp._raphanistrum    15 0.05025509  0.001   12
#Heracleum_sphondylium                        16 0.49177912  0.001    7
#Conopodium_majus                             17 0.58394048  0.001  182
#Anthriscus_sylvestris                        18 0.61245737  0.001  918
#Pteridium_aquilinum                          20 0.85154105  0.001   49
#Rosa_canina_agg.                             20 0.08151819  0.001   16
#Teucrium_scorodonia                          20 0.01972895  0.042  439
#Cirsium_arvense                              21 0.49335785  0.001  154
#Galium_aparine                               22 0.32095890  0.001  348
#Crataegus_monogyna                           23 0.65056506  0.001 1834
#Dioscorea_communis                           23 0.02726671  0.013    5
#Urtica_dioica                                24 0.36997392  0.001   18
#Equisetum_arvense                            25 0.96010826  0.001   17
#Prunella_vulgaris                            25 0.05567443  0.001    1
#Senecio_jacobaea                             25 0.02580840  0.037  221
#Epilobium_parviflorum                        25 0.02580836  0.021    4
#Festuca_ovina                                25 0.02486710  0.027    2
#Cynosurus_cristatus                          25 0.01414236  0.047    1


#To produce Indval results for the floristic tables in the supplementary
#information, the following code was used:

clustVEGCLUST <- cut(VEGCLUST$VEGCLUST_m1.1_dn0.90, 27, labels=FALSE)
VEGCLUSTIndval_1 <- indval(as.matrix(Vegclust_Plants[selection_criteria]), clustVEGCLUST + clustVEGCLUST, type= "long", allscores = TRUE)
summary(VEGCLUSTIndval_1, p = 1, too.many = 1000)

gr <- VEGCLUSTIndval_1$maxcls[VEGCLUSTIndval_1$pval<=1]
iv <- VEGCLUSTIndval_1$indcls[VEGCLUSTIndval_1$pval<=1]
pv <- VEGCLUSTIndval_1$pval[VEGCLUSTIndval_1$pval<=1]
fr <- apply(Vegclust_Plants[selection_criteria][,-1]>0, 2, sum)[VEGCLUSTIndval_1$pval<=1]
VEGCLUSTIndvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
VEGCLUSTIndvalsummary <- VEGCLUSTIndvalsummary[order(VEGCLUSTIndvalsummary$group, 
                                                     -VEGCLUSTIndvalsummary$indval),]
print(VEGCLUSTIndvalsummary)

write.csv(VEGCLUSTIndvalsummary, "FULL_INDVAL_RESULTS_25_VEGCLUST_COMMS.csv")

################################################################################
## 5. BIODIVERSITY INDEX CALCULATION
##
## Using hillR to generate biodiversity indices for the FSE data
################################################################################

#The next few lines of code perform a separate set of calculations, to investigate the
#mean, min, max and mode of different biodiversity indicators for each of the plant communities
#determined by previous Vegclust analysis.

#Preparing and checking the data...

Plant_Variables_Biodiversity_with_NA <- BMV_Vegclust %>% select(VEGCLUST_m1.1_dn0.90, q0, q1, q2)

Vegclust_Variables_Biodiversity <- Plant_Variables_Biodiversity_with_NA %>% drop_na(VEGCLUST_m1.1_dn0.90, q0, q1, q2)

head(Vegclust_Variables_Biodiversity)

#Calculating the mean, min, max and mode values of species richness firstly...

Vegclust_Variables_Biodiversity_Species_Richness_Mean <- Vegclust_Variables_Biodiversity %>% dplyr::group_by(VEGCLUST_m1.1_dn0.90)%>% summarise_at(vars(q0), list(name=mean)) 

Vegclust_Variables_Biodiversity_Species_Richness_Min <- Vegclust_Variables_Biodiversity %>% dplyr::group_by(VEGCLUST_m1.1_dn0.90)%>% summarise_at(vars(q0), list(name=min))   

Vegclust_Variables_Biodiversity_Species_Richness_Max <- Vegclust_Variables_Biodiversity %>% dplyr::group_by(VEGCLUST_m1.1_dn0.90)%>% summarise_at(vars(q0), list(name=max)) 

#This line of code produced the value "8", the mode.

names(sort(-table(Vegclust_Variables_Biodiversity$q0)))[1]

#Saving the data to put together in one file later...

write.csv(Vegclust_Variables_Biodiversity_Species_Richness_Mean, "Vegclust_Results_Dupes_Removed_Mean_Species_Richness.csv")

write.csv(Vegclust_Variables_Biodiversity_Species_Richness_Min, "Vegclust_Results_Dupes_Removed_Min_Species_Richness.csv")

write.csv(Vegclust_Variables_Biodiversity_Species_Richness_Max, "Vegclust_Results_Dupes_Removed_Max_Species_Richness.csv")

#Calculating the mean, min, max and mode values of Shannon Diversity secondly...

Vegclust_Variables_Biodiversity_Shannon_Mean <- Vegclust_Variables_Biodiversity %>% dplyr::group_by(VEGCLUST_m1.1_dn0.90)%>% summarise_at(vars(q1), list(name=mean)) 

Vegclust_Variables_Biodiversity_Shannon_Min <- Vegclust_Variables_Biodiversity %>% dplyr::group_by(VEGCLUST_m1.1_dn0.90)%>% summarise_at(vars(q1), list(name=min))   

Vegclust_Variables_Biodiversity_Shannon_Max <- Vegclust_Variables_Biodiversity %>% dplyr::group_by(VEGCLUST_m1.1_dn0.90)%>% summarise_at(vars(q1), list(name=max)) 

#Saving the data to put together in one file later...

write.csv(Vegclust_Variables_Biodiversity_Shannon_Mean, "Vegclust_Results_Dupes_Removed_Mean_Shannon.csv")

write.csv(Vegclust_Variables_Biodiversity_Shannon_Min, "Vegclust_Results_Dupes_Removed_Min_Shannon.csv")

write.csv(Vegclust_Variables_Biodiversity_Shannon_Max, "Vegclust_Results_Dupes_Removed_Max_Shannon.csv")

#Calculating the mean, min and max values of Simpson's Diversity secondly...

Vegclust_Variables_Biodiversity_Simpsons_Mean <- Vegclust_Variables_Biodiversity %>% dplyr::group_by(VEGCLUST_m1.1_dn0.90)%>% summarise_at(vars(q2), list(name=mean)) 

Vegclust_Variables_Biodiversity_Simpsons_Min <- Vegclust_Variables_Biodiversity %>% dplyr::group_by(VEGCLUST_m1.1_dn0.90)%>% summarise_at(vars(q2), list(name=min))   

Vegclust_Variables_Biodiversity_Simpsons_Max <- Vegclust_Variables_Biodiversity %>% dplyr::group_by(VEGCLUST_m1.1_dn0.90)%>% summarise_at(vars(q2), list(name=max)) 

#Saving the data to put together in one file later...

write.csv(Vegclust_Variables_Biodiversity_Simpsons_Mean, "Vegclust_Results_Dupes_Removed_Mean_Simpsons.csv")

write.csv(Vegclust_Variables_Biodiversity_Simpsons_Min, "Vegclust_Results_Dupes_Removed_Min_Simpsons.csv")

write.csv(Vegclust_Variables_Biodiversity_Simpsons_Max, "Vegclust_Results_Dupes_Removed_Max_Simpsons.csv")


################################################################################
## 6. VEGCLASS — APPLY FSE COMMUNITIES TO FARM4BIO DATA
##
## Align the Farm4Bio data with the chosen FSE clusters
## IMPORTANTLY - requires the data to have the same species recorded
################################################################################


## Read in the Farm4Bio data
#(This is a wide dataset with species as columns and sites as rows)
## Make sure this is the version that has the column header names
Farm4BioData   <- read.csv("Farm4Bio_All_Vegclust_Data_FSE_Name_Order_Metadata_Vo1.csv",
                           header = TRUE)

## Retain the species columns (as above for the FSE data)
Farm4BioData <- Farm4BioData %>% select(Chenopodium_album:Ranunculus_ficaria)

## Normalise (again same approach as used previously)
Farm4Biochord  <- decostand(Farm4BioData, "normalize")

## Apply the FSE vegclust communities to the Farm4Bio data
Farm4Bio.nc <- vegclass(FSE.nc.final, Farm4Biochord)


## Extract from the fitted object
groups_f4b <- defuzzify(Farm4Bio.nc)$cluster
table_f4b  <- table(groups_f4b)
print(table_f4b)

#Key result copied and pasted:
#groups_f4b
#F11 F12 F13 F14 F15 F16 F17 F19 F20 F22 F23 F25  F3  F4  F5  F8   N 
#100   2   1  31   4   1  20   8   5   7  27  39   1   8   3   3  77 

## Identify the assignment rate
n_f4b_assigned  <- sum(groups_f4b != "N")
pct_f4b_assigned <- round((n_f4b_assigned / nrow(Farm4Biochord)) * 100, 1)
cat("Farm4Bio quadrats assigned:", n_f4b_assigned, "/", nrow(Farm4Biochord),
    "(", pct_f4b_assigned, "%)\n")

#Key result copied and pasted:
#Farm4Bio quadrats assigned: 260 / 337 ( 77.2 %)

## Save as a file
write.csv(groups_f4b, file = "Farm4Bio_Vegclust_Groups.csv")

## Take a look at cluster variance
print(clustvar(Farm4Bio.nc))
print(clustvar(Farm4Bio.nc, defuzzify = TRUE))

#Hard results copied and pasted:
#M1        F2        F3        F4        F5        F6        F7        F8 
#NaN       NaN 0.6412511 0.5196357 0.6660480       NaN       NaN 0.7681248 
#F9       F10       F11       F12       F13       F14       F15       F16 
#NaN       NaN 0.6168594 0.6873744 0.5819471 0.5697241 0.6669084 0.6380022 
#F17       F18       F19       F20       F21       F22       F23       F24 
#0.6044452       NaN 0.5923266 0.7651407       NaN 0.6593427 0.5611943       NaN 
#F25 
#0.5468234 



#If you have had any issues with the code here, please contact Dr Savanna van Mesdag
#using either her Rothamsted email address, savanna.van-mesdag@rothamsted.ac.uk,
#or her current research profile, such as Researchgate.


################################################################################
## END
################################################################################


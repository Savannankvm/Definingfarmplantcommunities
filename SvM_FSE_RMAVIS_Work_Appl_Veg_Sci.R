#05 06 2026
#The code below was used to run RMAVIS analyses for the FSE data

#PLEASE BEAR IN MIND that the RShiny package RMAVIS has been updated since
#the code was originally run, so it may not be possible to exactly replicate
#the analyses using the exact code provided.

#install.packages("remotes")
install.packages(dplyr)

remotes::install_github("NERC-CEH/RMAVIS")

#Already installed.

library(remotes)
library(dplyr)

#Make sure to load remotes if the package needs to be installed.

#remotes::install_github("NERC-CEH/RMAVIS")

#This line of code installs the package, at least in theory, but there is another way of
#loading the package for use...

acceptedSpecies <- RMAVIS::acceptedSpecies

acceptedSpecies

write.csv(acceptedSpecies, "acceptedSpecies.csv")

?RMAVIS::similarityJaccard

?RMAVIS::speciesrichness

??RMAVIS

#Use this line of code to look up and check relevant documentation which will be needed
#for initial analyses.

RMAVIS::example_data[["Parsonage Down"]]

#This line of code runs an example, it might be useful to check this.

data <- RMAVIS::similarityJaccard(samp_df = RMAVIS::example_data[["Parsonage Down"]], 
                          comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                           habitatRestriction = c("CG"), 
                                                           col_name = "Pid3"),
                          samp_species_col = "Species", 
                          comp_species_col = "species",
                          samp_group_name = "Quadrat", 
                          comp_group_name = "Pid3",
                          comp_groupID_name = "NVC", 
                          remove_zero_matches = TRUE, 
                          average_comp = TRUE)

print(data)

dataParsonage <- RMAVIS::example_data[["Parsonage Down"]]

print(dataParsonage)

#Loading in the Beet data...
#This data file contains records for Beet fields specifically,
#as running all data would have taken much processing power for the computer
#used at the time. With smaller datasets or with more powerful computers
#running the entire dataset is more feasible.
#Species names were standardised according to accepted species RMAVIS list, 
#with one column containing species and each row containing individual species 
#records within a quadrat sample.


BEETdata_corrected_raw <- read.csv("BMV_Data_Long_Format_RMAVIS_BEET_Corrected_Vo4.csv",  
                         header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(BEETdata_corrected_raw)

#This is the line of code for the NVC analyses. Note that the "c("CG") part of the code makes
#this analysis specific to CG habitats. It will be necessary to include other habitats
#such as OV etc that are most likely to be found in arable field margins, partially
#based on the communities that were identified in the preliminary analyses.
#Also make sure to change the "sampe_species_col" column to your "Species" field and the 
#"sample_group_name" column to your quadrat column/similar.


BEETdata_corrected <- RMAVIS::similarityJaccard(samp_df = BEETdata_corrected_raw, 
                                      comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                       habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                       col_name = "Pid3"),
                                      samp_species_col = "SPECIES", 
                                      comp_species_col = "species",
                                      samp_group_name = "QUADRAT", 
                                      comp_group_name = "Pid3",
                                      comp_groupID_name = "NVC", 
                                      remove_zero_matches = TRUE, 
                                      average_comp = TRUE)


write.csv(BEETdata_corrected, "BEET.data_corrected_portion.csv")



library(dplyr)

BEETdata_corrected <- read.csv("BEET.data_corrected_portion.csv",  
                     header = TRUE)

B2 <- BEETdata_corrected %>%
  group_by(QUADRAT) %>%
  summarise(max = max(Similarity, na.rm=TRUE))

print(B2, n=15)

options(dplyr.print_max = 1e9)

print(B2, n = 51)


print(B1, n = 500)


#Loading in the Maize data...
#This data file contains records for Maize fields specifically,
#as running all data would have taken much processing power for the computer
#used at the time. With smaller datasets or with more powerful computers
#running the entire dataset is more feasible.
#Species names were standardised according to accepted species RMAVIS list, 
#with one column containing species and each row containing individual species 
#records within a quadrat sample.


MAIZEdata_raw <- read.csv("BMV_RMAVIS_MAIZE_Vo1.csv",  
                         header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(MAIZEdata_raw)



MAIZEdata <- RMAVIS::similarityJaccard(samp_df = MAIZEdata_raw, 
                                      comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                       habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                       col_name = "Pid3"),
                                      samp_species_col = "SPECIES", 
                                      comp_species_col = "species",
                                      samp_group_name = "QUADRAT", 
                                      comp_group_name = "Pid3",
                                      comp_groupID_name = "NVC", 
                                      remove_zero_matches = TRUE, 
                                      average_comp = TRUE)


write.csv(MAIZEdata, "MAIZEdata.csv")


library(dplyr)

MAIZEdata <- read.csv("MAIZEdata.csv",  
                               header = TRUE)

M1 <- MAIZEdata %>%
  group_by(QUADRAT) %>%
  summarise(max = max(Similarity, na.rm=TRUE))

print(M1, n=15)

options(dplyr.print_max = 1e9)

print(M1, n = 500)

#12th June, latest one, 1002313.

print (M1)


#Loading in the Spring Oil Seed Rape data...
#This data file contains records for Spring OSR fields specifically,
#as running all data would have taken much processing power for the computer
#used at the time. With smaller datasets or with more powerful computers
#running the entire dataset is more feasible.
#Species names were standardised according to accepted species RMAVIS list, 
#with one column containing species and each row containing individual species 
#records within a quadrat sample.

SPRINGOSRdata_raw <- read.csv("BMV_SPRING_OILSEED_RAPE_RMAVIS_Data_1.csv",  
                          header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(SPRINGOSRdata_raw)

#Running the analysis and generating the output needed, with similarity values...

SPRINGOSRdata <- RMAVIS::similarityJaccard(samp_df = SPRINGOSRdata_raw, 
                                       comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                        habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                        col_name = "Pid3"),
                                       samp_species_col = "SPECIES", 
                                       comp_species_col = "species",
                                       samp_group_name = "QUADRAT", 
                                       comp_group_name = "Pid3",
                                       comp_groupID_name = "NVC", 
                                       remove_zero_matches = TRUE, 
                                       average_comp = TRUE)


#Helen helped me with this - instead of deleting values manually, I have now filtered
#all of the highest similarity values for each quadrat, for the spring oilseed rape data

write.csv(SPRINGOSRdata, "SPRINGOSRdata.csv")

SPRINGOSRdata <- read.csv("SPRINGOSRdata.csv",  
                      header = TRUE)

SPOSR1 <- SPRINGOSRdata %>%
  group_by(QUADRAT) %>%
  filter(Similarity == max(Similarity))

print(SPOSR1)

write.csv(SPOSR1, "SPRINGOSRdataVo1.csv")

#Loading in the Winter Oil Seed Rape data...
#This data file contains records for Winter OSR fields specifically,
#as running all data would have taken much processing power for the computer
#used at the time. With smaller datasets or with more powerful computers
#running the entire dataset is more feasible.
#Species names were standardised according to accepted species RMAVIS list, 
#with one column containing species and each row containing individual species 
#records within a quadrat sample.

WINTEROSRdata <- read.csv("WINTEROSRdata.csv",  
                          header = TRUE)

WOSR1 <- WINTEROSRdata %>%
  group_by(QUADRAT) %>%
  filter(Similarity == max(Similarity))

print(WOSR1)

write.csv(WOSR1, "WINTEROSRdataVo1.csv")


#Now for data for winter oilseed rape...

#Loading the raw data

WINTEROSRdata_raw <- read.csv("BMV_WINTER_OILSEED_RAPE_RMAVIS_Data_1.csv",  
                              header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(WINTEROSRdata_raw)

#Running the analysis and generating the output needed, with similarity values...

WINTEROSRdata <- RMAVIS::similarityJaccard(samp_df = WINTEROSRdata_raw, 
                                       comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                        habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                        col_name = "Pid3"),
                                       samp_species_col = "SPECIES", 
                                       comp_species_col = "species",
                                       samp_group_name = "QUADRAT", 
                                       comp_group_name = "Pid3",
                                       comp_groupID_name = "NVC", 
                                       remove_zero_matches = TRUE, 
                                       average_comp = TRUE)


write.csv(WINTEROSRdata, "WINTEROSRdata.csv")

###############################

#####Rerunning the analyses, but removing woody species######


BEETdata_raw_not_woody <- read.csv("BMV_Data_Long_Format_RMAVIS_BEET_WITHOUT_WOODY_INPUT_Vo4.csv",  
                         header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(BEETdata_raw_not_woody)


BEETdata_not_woody <- RMAVIS::similarityJaccard(samp_df = BEETdata_raw_not_woody, 
                                      comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                       habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                       col_name = "Pid3"),
                                      samp_species_col = "SPECIES", 
                                      comp_species_col = "species",
                                      samp_group_name = "QUADRAT", 
                                      comp_group_name = "Pid3",
                                      comp_groupID_name = "NVC", 
                                      remove_zero_matches = TRUE, 
                                      average_comp = TRUE)

write.csv(BEETdata_not_woody, "FSE_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv")

BEETdata_not_woody <- read.csv("FSE_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv",  
                          header = TRUE)

library(dplyr)

BEET_NOT_WOODY1 <- BEETdata_not_woody %>%
  group_by(QUADRAT) %>%
  filter(Similarity == max(Similarity))

print(BEET_NOT_WOODY1)

write.csv(BEET_NOT_WOODY1, "FSE_RMAVIS_Results_Long_Format_OUTPUT_Max_Similarities.csv")



#Now to do the same for the Maize data...


MAIZEdata_raw_not_woody <- read.csv("FSE_MAIZE_RMAVIS_Long_Format_WITHOUT_WOODY_INPUT_Vo1.csv",  
                                   header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(MAIZEdata_raw_not_woody)


MAIZEdata_not_woody <- RMAVIS::similarityJaccard(samp_df = MAIZEdata_raw_not_woody, 
                                                comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                                 habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                                 col_name = "Pid3"),
                                                samp_species_col = "SPECIES", 
                                                comp_species_col = "species",
                                                samp_group_name = "QUADRAT", 
                                                comp_group_name = "Pid3",
                                                comp_groupID_name = "NVC", 
                                                remove_zero_matches = TRUE, 
                                                average_comp = TRUE)

write.csv(MAIZEdata_not_woody, "FSE_MAIZE_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv")

MAIZEdata_not_woody <- read.csv("FSE_MAIZE_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv",  
                               header = TRUE)

library(dplyr)

MAIZE_NOT_WOODY1 <- MAIZEdata_not_woody %>%
  group_by(QUADRAT) %>%
  filter(Similarity == max(Similarity))

print(MAIZE_NOT_WOODY1)

write.csv(MAIZE_NOT_WOODY1, "FSE_MAIZE_RMAVIS_Results_Long_Format_OUTPUT_Max_Similarities.csv")

#Now to run the Spring OSR data...

SPRINGOSRdata_raw_not_woody <- read.csv("FSE_SPRING_OSR_RMAVIS_Long_Format_INPUT_WITHOUT_WOODY_Vo1.csv",  
                                    header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(SPRINGOSRdata_raw_not_woody)


SPRINGOSRdata_not_woody <- RMAVIS::similarityJaccard(samp_df = SPRINGOSRdata_raw_not_woody, 
                                                 comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                                  habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                                  col_name = "Pid3"),
                                                 samp_species_col = "SPECIES", 
                                                 comp_species_col = "species",
                                                 samp_group_name = "QUADRAT", 
                                                 comp_group_name = "Pid3",
                                                 comp_groupID_name = "NVC", 
                                                 remove_zero_matches = TRUE, 
                                                 average_comp = TRUE)

write.csv(SPRINGOSRdata_not_woody, "FSE_SPRINGOSR_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv")

SPRINGOSRdata_not_woody <- read.csv("FSE_SPRINGOSR_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv",  
                                header = TRUE)

library(dplyr)

SPRINGOSR_NOT_WOODY1 <- SPRINGOSRdata_not_woody %>%
  group_by(QUADRAT) %>%
  filter(Similarity == max(Similarity))

print(SPRINGOSR_NOT_WOODY1)

write.csv(SPRINGOSR_NOT_WOODY1, "FSE_SPRINGOSR_RMAVIS_Results_Long_Format_OUTPUT_Max_Similarities.csv")

#Now to run the Winter OSR data...

WINTEROSRdata_raw_not_woody <- read.csv("FSE_WINTER_OSR_RMAVIS_Long_Format_INPUT_WITHOUT_WOODY_Vo1.csv",  
                                        header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(WINTEROSRdata_raw_not_woody)


WINTEROSRdata_not_woody <- RMAVIS::similarityJaccard(samp_df = WINTEROSRdata_raw_not_woody, 
                                                     comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                                      habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                                      col_name = "Pid3"),
                                                     samp_species_col = "SPECIES", 
                                                     comp_species_col = "species",
                                                     samp_group_name = "QUADRAT", 
                                                     comp_group_name = "Pid3",
                                                     comp_groupID_name = "NVC", 
                                                     remove_zero_matches = TRUE, 
                                                     average_comp = TRUE)

write.csv(WINTEROSRdata_not_woody, "FSE_WINTEROSR_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv")

WINTEROSRdata_not_woody <- read.csv("FSE_WINTEROSR_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv",  
                                    header = TRUE)

library(dplyr)

WINTEROSR_NOT_WOODY1 <- WINTEROSRdata_not_woody %>%
  group_by(QUADRAT) %>%
  filter(Similarity == max(Similarity))

print(WINTEROSR_NOT_WOODY1)

write.csv(WINTEROSR_NOT_WOODY1, "FSE_WINTEROSR_RMAVIS_Results_Long_Format_OUTPUT_Max_Similarities.csv")

#I'm going to redo the Beet analysis for the woody and non-woody species combined:

BEETdata_raw <- read.csv("BMV_Data_Long_Format_RMAVIS_BEET_INPUT_Vo3.csv",  
                                   header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(BEETdata_raw)


BEETdata <- RMAVIS::similarityJaccard(samp_df = BEETdata_raw, 
                                                comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                                 habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                                 col_name = "Pid3"),
                                                samp_species_col = "SPECIES", 
                                                comp_species_col = "species",
                                                samp_group_name = "QUADRAT", 
                                                comp_group_name = "Pid3",
                                                comp_groupID_name = "NVC", 
                                                remove_zero_matches = TRUE, 
                                                average_comp = TRUE)

write.csv(BEETdata, "FSE_Beet_WITH_WOODY_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv")

BEETdata <- read.csv("FSE_Beet_WITH_WOODY_RMAVIS_Results_Long_Format_OUTPUT_Vo1.csv",  
                               header = TRUE)

library(dplyr)

BEET1 <- BEETdata %>%
  group_by(QUADRAT) %>%
  filter(Similarity == max(Similarity))

print(BEET1)

write.csv(BEET1, "FSE_Beet_WITH_WOODY_RMAVIS_Results_Long_Format_OUTPUT_Max_Similarities.csv")





#RMAVIS::master_data

#You probably won't need to run this again - but this is where the Ellenberg data
#come from



#savannahs_data |> dplyr::left_join
#(RMAVIS::master_data, by = "species") |> dplyr::group_by("transect_name") |> dplyr::summarise
#("F_mean" = mean(F, na.rm = TRUE, "L_mean" = mean(L, na.rm = TRUE)) |> dplyr::ungroup()

#This line of code/a very similar line of code would provide Ellenberg values for the different
#species in the communities.

#Loading the community attribute data...

ca <- RMAVIS:: nvc_community_attributes

write.csv(ca, "NVC_Community_Attribute_Table_RMAVIS.csv")

#I now need to generate species richness values for each without-woody quadrat...


BEETdata_raw_not_woody <- read.csv("BMV_Data_Long_Format_RMAVIS_BEET_WITHOUT_WOODY_INPUT_Vo4.csv",  
                                   header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(BEETdata_raw_not_woody)

library(dplyr)

BEETdata_raw_not_woody_species_richness <- BEETdata_raw_not_woody %>% dplyr::group_by(QUADRAT)%>% summarise(n=dplyr::n_distinct(SPECIES))

write.csv(BEETdata_raw_not_woody_species_richness, "BMV_Data_Long_Format_RMAVIS_BEET_WITHOUT_WOODY_INPUT_Species_richness_Vo1.csv")

MAIZEdata_raw_not_woody <- read.csv("FSE_MAIZE_RMAVIS_Long_Format_WITHOUT_WOODY_INPUT_Vo1.csv",  
                                    header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(MAIZEdata_raw_not_woody)

MAIZEdata_raw_not_woody_species_richness <- MAIZEdata_raw_not_woody %>% dplyr::group_by(QUADRAT)%>% summarise(n=dplyr::n_distinct(SPECIES))

write.csv(MAIZEdata_raw_not_woody_species_richness, "FSE_MAIZE_RMAVIS_Long_Format_WITHOUT_WOODY_INPUT_Species_richness_Vo1.csv")

SPRINGOSRdata_raw_not_woody <- read.csv("FSE_SPRING_OSR_RMAVIS_Long_Format_INPUT_WITHOUT_WOODY_Vo1.csv",  
                                        header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(SPRINGOSRdata_raw_not_woody)

SPRINGOSRdata_raw_not_woody_species_richness <- SPRINGOSRdata_raw_not_woody %>% dplyr::group_by(QUADRAT)%>% summarise(n=dplyr::n_distinct(SPECIES))

write.csv(SPRINGOSRdata_raw_not_woody_species_richness, "FSE_SPRING_OSR_RMAVIS_Long_Format_INPUT_WITHOUT_WOODY_Species_richness_Vo1.csv")

WINTEROSRdata_raw_not_woody <- read.csv("FSE_WINTER_OSR_RMAVIS_Long_Format_INPUT_WITHOUT_WOODY_Vo1.csv",  
                                        header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(WINTEROSRdata_raw_not_woody)

WINTEROSRdata_raw_not_woody_species_richness <- WINTEROSRdata_raw_not_woody %>% dplyr::group_by(QUADRAT)%>% summarise(n=dplyr::n_distinct(SPECIES))

write.csv(WINTEROSRdata_raw_not_woody_species_richness, "FSE_WINTER_OSR_RMAVIS_Long_Format_INPUT_WITHOUT_WOODY_Species_richness_Vo1.csv")


#Loading in latest version of full data



BMV_All <- read.csv("FSE_Combined_Data_Vo7.csv",  
                                    header = TRUE)



BMV_All_Species_Richness_Mean <- BMV_All %>% dplyr::group_by(NVC_without_woody)%>% summarise_at(vars(Without_Woody_Species_Richness), list(name=mean))   

BMV_All_Species_Richness_Min <- BMV_All %>% dplyr::group_by(NVC_without_woody)%>% summarise_at(vars(Without_Woody_Species_Richness), list(name=min)) 

BMV_All_Species_Richness_Max <- BMV_All %>% dplyr::group_by(NVC_without_woody)%>% summarise_at(vars(Without_Woody_Species_Richness), list(name=max)) 

BMV_All_Species_Richness_Sum <- BMV_All %>% dplyr::group_by(NVC_without_woody)%>% summarise_at(vars(Without_Woody_Species_Richness), list(name=sum)) 

BMV_All_Species_Richness_Mode <- BMV_All %>% dplyr::group_by(NVC_without_woody)%>% summarise_at(vars(Without_Woody_Species_Richness), list(name=mode)) 


write.csv(BMV_All_Species_Richness_Mean, "RMAVIS_Results_Mean_Species_Richness_Without_Woody.csv")

write.csv(BMV_All_Species_Richness_Min, "RMAVIS_Results_Min_Species_Richness_Without_Woody.csv")

write.csv(BMV_All_Species_Richness_Max, "RMAVIS_Results_Max_Species_Richness_Without_Woody.csv")

write.csv(BMV_All_Species_Richness_Sum, "RMAVIS_Results_Sum_Species_Richness_Without_Woody.csv")

write.csv(BMV_All_Species_Richness_Mode, "RMAVIS_Results_Mode_Species_Richness_Without_Woody.csv")


get_mode <- function(x) {
  
  unique_values <- unique(x)
  
  frequencies <- tabulate(match(x, unique_values))
  
  unique_values[which.max(frequencies)]
  
}

BMV_All %>% group_by(NVC_without_woody) %>% summarise(mode = get_mode(Without_Woody_Species_Richness))




#05/06/2026

#S van Mesdag and Z Marshall

#The following lines of code were used to run RMAVIS to generate NVC results
#for the Farm4Bio data

#PLEASE BEAR IN MIND that the RShiny package has been updated since this code has
#been run so it may not be possible to run this code without adpating things
#to reflect recent changes in the RShiny package.

install.packages("remotes")

library(remotes)

#Make sure to load remotes if the package needs to be installed.

remotes::install_github("NERC-CEH/RMAVIS")

#This line of code installs the package, at least in theory, but there is another way of
#loading the package for use...

acceptedSpecies <- RMAVIS::acceptedSpecies

subset_vcData <- RMAVIS::subset_vcData

acceptedSpecies

write.csv(acceptedSpecies, "acceptedSpecies.csv")

?RMAVIS::similarityJaccard

#Use this line of code to look up and check relevant documentation which will be needed
#for initial analyses.

RMAVIS::example_data[["Parsonage Down"]]

#This line of code runs an example, it might be useful to check this.

#In this particular example, I now need to pivot_longer the plant data in the table, so all of the
#species fit into one column labelled "SPECIES" and all of the cover values fit into a column
#labelled "COVER."

#install.packages("devtools")

library(devtools)

#install.packages("permute")

library(permute)

library(tidyr)

library(snakecaser)

library(dplyr)

#Loading the data....
#This data file contained species named standardised according to
#accepted species RMAVIS list, with one column containing species
#and each row containing individual species records within a quadrat sample.

FARM4BIOdata_raw <- read.csv("Farm4Bio_RMAVIS_Data_1.csv",  
                         header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

help(pivot_longer)

#Attempting pivot_longer...

FARM4BIOdata_raw %>% pivot_longer(cols = (Acer campestre:Viola arvensis), 
                                          names_to = c("SPECIES", "COVER"),
                                          values_to = "count")

FARM4BIOdata_long = FARM4BIOdata_raw %>% pivot_longer(cols = (Acer.campestre:Viola.arvensis), 
                                                      names_to = c("SPECIES"),
                                                      values_to = "count")

FARM4BIOdata_long <- write.csv()


write.csv(FARM4BIOdata_long, "FARM4BIOdata_long.csv")

#I have edited the file within Excel to make sure the names are correct for the analysis
#and to remove all species where the cover values = 0.

FARM4BIOdata_raw <- read.csv("FARM4BIOdata_long.csv",  
                             header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

print(FARM4BIOdata_raw)


FARM4BIOdata_raw <- RMAVIS::similarityJaccard(samp_df = FARM4BIOdata_raw, 
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

write.csv(FARM4BIOdata_raw, "FARM4BIOdata_RMAVIS.csv")

#Filtering values so only the highest for each quadrat remain

library(dplyr)


FARM4BIOdata_raw <- read.csv("FARM4BIOdata_RMAVIS.csv",  
                          header = TRUE)

FARM4BIO1 <- FARM4BIOdata_raw %>%
  group_by(QUADRAT) %>%
  filter(Similarity == max(Similarity))

print(FARM4BIO1)

write.csv(FARM4BIO1, "Farm4BioRMAVISResultsVo1.csv")



#Running the same analysis and calculations, but with the 'woody'
#data removed, following the list of species set in the file titled:
#"RMAVIS_Species_Lists_Woody_and_not_woody.xlsx".
#####Please note that I am having trouble with this line of code######


FARM4BIOdata_WW <- read.csv("FARM4BIOdata_long_INPUT_RMAVIS_Without_Woody.csv",  
                             header = TRUE)|>
  
  tibble::as_tibble() |>
  
  dplyr::mutate("QUADRAT" = as.character(QUADRAT))

# Check which species are not accepted
setdiff(unique(FARM4BIOdata_WW$SPECIES), RMAVIS::accepted_taxa$taxon_name) # <--- 13 species not accepted by RMAVIS

# Check whether any cover estimates are supplied
isTRUE(all(!is.na(unique(FARM4BIOdata_WW$COVER))))

# Check whether all cover estimates are supplied
isTRUE(all(!is.na(FARM4BIOdata_WW$COVER)))

# Check whether there is any missing data in the Year column
isTRUE(all(!is.na(FARM4BIOdata_WW$DATE)))

# Check whether there is any missing data in the Group column
isTRUE(all(FARM4BIOdata_WW$GROUP != ""))

# Fix data ----------------------------------------------------------------

# Here I adjust the species not found in the RMAVIS::accepted_taxa using info 
# UKVegTB::taxa_lookup, which provides a lookup between the taxa accepted by
# RMAVIS and their synonyms as present in the UKSI.
FARM4BIOdata_WW_fixed <- FARM4BIOdata_WW |>
  dplyr::mutate(
    "SPECIES" = dplyr::case_when(
      SPECIES == "Picris echioides" ~ "Helminthotheca echioides",
      SPECIES == "Taraxacum officinale agg." ~ "Taraxacum agg.",
      SPECIES == "Bromus commutatus" ~ "Bromus racemosus subsp. commutatus",
      SPECIES == "Bromus sterilis" ~ "Anisantha sterilis",
      SPECIES == "Elytrigia repens" ~ "Elymus repens",
      
      SPECIES == "Anagallis arvensis" ~ "Lysimachia arvensis",
      SPECIES == "Polygonum aviculare" ~ "Polygonum aviculare sensu Stace",
      SPECIES == "Centaurea nigra sensu stricto" ~ "Centaurea nigra",
      SPECIES == "Hypericum perforatum" ~ "Hypericum perforatum sensu Stace",
      SPECIES == "Chamerion angustifolium" ~ "Chamaenerion angustifolium",
      
      SPECIES == "Galium mollugo" ~ "Galium album",
      SPECIES == "Vicia tetrasperma" ~ "Ervum tetraspermum",
      SPECIES == "Senecio jacobaea subsp. jacobaea" ~ "Jacobaea vulgaris subsp. vulgaris",
      
      TRUE ~ SPECIES
    )
  )

# Check which species are not accepted, again
setdiff(unique(FARM4BIOdata_WW_fixed$SPECIES), RMAVIS::accepted_taxa$taxon_name) # <--- 0 species not accepted by RMAVIS

#All OK, all of the relevant species have been changed appropriately.

# Calculate the Jaccard similarities --------------------------------------
dat_jac_sims <- RMAVIS::similarityJaccard(samp_df = FARM4BIOdata_WW_fixed, 
                                          comp_df = RMAVIS::subset_nvcData(nvc_data = RMAVIS::nvc_pquads, 
                                                                          habitatRestriction = c("CG", "H", "MG", "A", "SD", "OV"), 
                                                                          col_name = "psq_id"),
                                          samp_species_col = "SPECIES", 
                                          comp_species_col = "nvc_taxon_name",
                                          samp_group_name = "QUADRAT", 
                                          comp_group_name = "psq_id",
                                          comp_groupID_name = "nvc_code", 
                                          remove_zero_matches = TRUE, 
                                          average_comp = TRUE) |>
  tibble::as_tibble()



write.csv(dat_jac_sims, "FARM4BIOdata_WW_RMAVIS.csv")

#Filtering values so only the highest for each quadrat remain

library(dplyr)


FARM4BIOdata_WW_raw <- read.csv("FARM4BIOdata_WW_RMAVIS.csv",  
                             header = TRUE)

FARM4BIO1 <- FARM4BIOdata_WW_raw %>%
  group_by(QUADRAT) %>%
  filter(Similarity == max(Similarity))

print(FARM4BIO1)

write.csv(FARM4BIO1, "Farm4BioRMAVISResultsWWVo1.csv")




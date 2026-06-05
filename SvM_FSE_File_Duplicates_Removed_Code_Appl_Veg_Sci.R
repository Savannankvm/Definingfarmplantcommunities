#05 06 2026
#S van Mesdag

#This code was used to identify and remove duplicates from the original FSE data.
#The data file was then incorporated into the final METADATA file in the repository.

#Set directory as needed

#Loading in the data
#This file has species as columns and each row represents a quadrat sample

FSEDataWide <- read.csv("BMV_Selected_Metadata_Wide_Vo1.csv",  
                        header = TRUE)

glimpse(FSEDataWide)

str(FSEDataWide)

unique(FSEDataWide$Species) |> sort()
unique(FSEDataWide$QUADRAT) |> length()
unique(FSEDataWide$FEATURE_CODE) |> sort()

summary(FSEDataWide$Cover)

## Identify duplicates ----------------------------------------------------
raw_dat_dupes <- FSEDataWide |>
  dplyr::group_by(QUADRAT, VISIT_DATE) |>
  dplyr::filter(dplyr::n() > 1) |>
  dplyr::ungroup() |>
  dplyr::select(QUADRAT, VISIT_DATE) |>
  dplyr::arrange(QUADRAT)

print(raw_dat_dupes)

#Creating new file with duplicates removed

FSEDataWideUnique <- FSEDataWide[!duplicated(FSEDataWide$QUADRAT),]

write.csv (FSEDataWideUnique, file = "FSEDataWideUnique_Vo1.csv")

#If you have any issues with running this code, please contact Dr Savanna van Mesdag
#using her Rothamsted email address, savanna.van-mesdag@rothamsted.ac.uk, or 
#her current Research profile, such as Researchgate.

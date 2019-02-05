rm(list = ls())

setwd("~/Documents/hotterbetterprokaryotes/Code/R_Scripts/")

pathogen_data <- read.csv("../../Data/SpeciesInteractions_EID2.csv")
prokaryote_data <- read.csv("../../Data/summaries/non_aggregate_data_fix.csv")

# add pathogen columns to the summary
prokaryote_data$Pathogen <- NA
prokaryote_data$Pathogen_type <- NA

# I want to check whether species in the database are also in the pathogen database and if so, assign them the pathogen label

for (i in 1:length(prokaryote_data$Species)){
  # lets paste genus and species together, so we can avoid strain names
  species <- paste(as.character(prokaryote_data$ConGenus[i]), as.character(prokaryote_data$ConSpecies[i]), sep = " ")
  is_pathogen <- grep(species, pathogen_data$Cargo, ignore.case=TRUE)

  if (length(is_pathogen) > 0){
    prokaryote_data$Pathogen[i] <- "Pathogen"
    prokaryote_data$Pathogen_type[i] <- unique(as.character(pathogen_data[is_pathogen,]$Carrier.classification))[1]
  }
  else{
    prokaryote_data$Pathogen[i] <- "Non-pathogen"
  }
}

# and add in an extra catch where bacteria supposedly infecting bacteria are not listed as pathogens:
prokaryote_data[prokaryote_data$Pathogen_type == "Bacteria" & !is.na(prokaryote_data$Pathogen_type),]$Pathogen <- "Non-pathogen"
prokaryote_data[prokaryote_data$Pathogen_type == "Bacteria" & !is.na(prokaryote_data$Pathogen_type),]$Pathogen_type <- NA

write.csv(prokaryote_data, "../../Data/summaries/summary_pathogens.csv", row.names = FALSE)

# are pathogens different to others?

mean_ES <- mean(prokaryote_data$E, na.rm=TRUE)

mean_pathogens <- mean(subset(prokaryote_data, Pathogen == "Pathogen")$E, na.rm=TRUE)
mean_non_pathogens <- mean(subset(prokaryote_data, Pathogen == "Non-pathogen")$E, na.rm=TRUE)

good_data <- prokaryote_data[prokaryote_data$Trait == "Specific Growth Rate" & prokaryote_data$Points_Before_Peak > 2,]
ggplot(good_data, aes(x = E)) + geom_density(aes(fill = Pathogen), alpha = 0.7) + xlim(0, 5)

ggplot(good_data[good_data$Pathogen == "Pathogen",], aes(x = E)) + geom_density(aes(fill = Pathogen_type), alpha = 0.7) + xlim(0, 2)

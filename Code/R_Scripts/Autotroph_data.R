## Check taxonomy of Sofia's data by comparing back to BioTraits

rm(list = ls())

setwd("~/Documents/hotterbetterprokaryotes/Code/R_Scripts/")

# first combine Sofia's aquatic and terrestrial results

# aq_data <- read.csv("../../Data/sofia_data/resultsAQ.csv")
# ter_data <- read.csv("../../Data/sofia_data/resultsTER.csv")
# 
# aq_data$habitat <- "Aquatic"
# ter_data$habitat <- "Terrestrial"
# 
# full_data <- rbind(aq_data, ter_data)
# 
# # only interested in respiration
# autotroph_data <- full_data[full_data$id_process == "respiration rate",]
# autotroph_data <- autotroph_data[!is.na(autotroph_data$E_sch),]
# 
# # also only want fits with sensible E
# autotroph_data <- autotroph_data[autotroph_data$E_sch <= 3,]
# 
# # load biotraits data
# 
# biotraits <- read.csv("~/Documents/biotraitsdb/data/GlobalDataset_v0.5.csv")
# 
# autotroph_ids <- as.character(unique(autotroph_data$id))
# 
# autotroph_biotraits <- biotraits[biotraits$originalid %in% autotroph_ids,]
# 
# length(unique(autotroph_biotraits$originalid)) # length checks out, we're good
# 
# unique(autotroph_biotraits$interactor1kingdom)
# 
# # remove (cyano)bacteria so we arent counting them twice!!
# 
# autotroph_biotraits <- autotroph_biotraits[autotroph_biotraits$interactor1kingdom != "Bacteria",]
# 
# # could if we wanted write this dataset out and do the whole fitting ourselves with the python code, but not necessary for now
# # now match these IDs back to the original data to remove bacteria
# 
# non_bact_IDs <- as.character(unique(autotroph_biotraits$originalid))
# 
# autotroph_data <- autotroph_data[autotroph_data$id %in% non_bact_IDs,]
# 
# # write this out for use elsewhere
# #write.csv(autotroph_data, "../../Data/sofia_data/DataAutotrophs.csv", row.names = FALSE)
# 
# ####### now back to checking whats actually in sofias data! #######
# 
# length(unique(autotroph_biotraits$interactor1species))
# 
# length(unique(autotroph_biotraits[autotroph_biotraits$interactor1kingdom == "Chromista",]$interactor1species))
# unique(autotroph_biotraits[autotroph_biotraits$interactor1kingdom == "Chromista",]$interactor1class) # brown algae
# 
# unique(autotroph_biotraits[autotroph_biotraits$interactor1kingdom == "Plantae",]$interactor1phylum)
# 
# # red algae
# length(unique(autotroph_biotraits[autotroph_biotraits$interactor1phylum == "Rhodophyta",]$interactor1species))
# 
# unique(autotroph_biotraits[autotroph_biotraits$interactor1phylum == "Streptophyta",]$interactor1class)
# 
# # green algae
# length(unique(autotroph_biotraits[autotroph_biotraits$interactor1phylum == "Chlorophyta",]$interactor1species))
# length(unique(autotroph_biotraits[autotroph_biotraits$interactor1class %in% c("Klebsormidiophyceae", "Zygnemophyceae"),]$interactor1species))
# 
# # mosses
# length(unique(autotroph_biotraits[autotroph_biotraits$interactor1class %in% c("Sphagnopsida", "Bryopsida"),]$interactor1species))
# 
# # vascular plants
# length(unique(autotroph_biotraits[autotroph_biotraits$interactor1class %in% c("Liliopsida", "Magnoliopsida", "Pinopsida"),]$interactor1species))
# 
# 
# #########
# # lets add taxonomy to the autotroph fitting results
# #
# 
# autotroph_data$ConKingdom <- NA
# autotroph_data$ConPhylum <- NA
# autotroph_data$ConClass <- NA
# autotroph_data$ConOrder <- NA
# autotroph_data$ConFamily <- NA
# autotroph_data$ConGenus <- NA
# autotroph_data$ConSpecies <- NA
# 
# for(i in 1:length(non_bact_IDs)){
#   autotroph_data$ConKingdom[i] <- 
#     as.character(unique(autotroph_biotraits[autotroph_biotraits$originalid == as.character(autotroph_data[i,]$id),]$interactor1kingdom))
#   autotroph_data$ConPhylum[i] <- 
#     as.character(unique(autotroph_biotraits[autotroph_biotraits$originalid == as.character(autotroph_data[i,]$id),]$interactor1phylum))
#   autotroph_data$ConClass[i] <- 
#     as.character(unique(autotroph_biotraits[autotroph_biotraits$originalid == as.character(autotroph_data[i,]$id),]$interactor1class))
#   autotroph_data$ConOrder[i] <- 
#     as.character(unique(autotroph_biotraits[autotroph_biotraits$originalid == as.character(autotroph_data[i,]$id),]$interactor1order))
#   autotroph_data$ConFamily[i] <- 
#     as.character(unique(autotroph_biotraits[autotroph_biotraits$originalid == as.character(autotroph_data[i,]$id),]$interactor1family))
#   autotroph_data$ConGenus[i] <- 
#     as.character(unique(autotroph_biotraits[autotroph_biotraits$originalid == as.character(autotroph_data[i,]$id),]$interactor1genus))
#   autotroph_data$ConSpecies[i] <- 
#     as.character(unique(autotroph_biotraits[autotroph_biotraits$originalid == as.character(autotroph_data[i,]$id),]$interactor1species))
# }
# 
# # and also add the "common name" higher level groups
# 
# autotroph_data$TaxonomicGroup <- NA
# 
# autotroph_data[autotroph_data$ConClass == "Phaeophyceae",]$TaxonomicGroup <- "Brown Algae"
# autotroph_data[autotroph_data$ConPhylum == "Rhodophyta",]$TaxonomicGroup <- "Red Algae"
# autotroph_data[autotroph_data$ConPhylum == "Chlorophyta",]$TaxonomicGroup <- "Green Algae"
# autotroph_data[autotroph_data$ConClass %in% c("Klebsormidiophyceae", "Zygnemophyceae"),]$TaxonomicGroup <- "Green Algae"
# autotroph_data[autotroph_data$ConClass %in% c("Sphagnopsida", "Bryopsida"),]$TaxonomicGroup <- "Mosses"
# autotroph_data[autotroph_data$ConClass %in% c("Liliopsida", "Magnoliopsida", "Pinopsida"),]$TaxonomicGroup <- "Vascular Plants"
# 
# # add an extra column sort of for the python stuff later
# autotroph_data$Trophy <- "Autotroph"
# 
# # now we should align the column names to our prokaryote data so it will go through the python code with no fuss
# names(autotroph_data)[5] <- "E"
# names(autotroph_data)[6] <- "E_std"
# 
# # now write this version out instead of the old one
# write.csv(autotroph_data, "../../Data/sofia_data/DataAutotrophs.csv", row.names = FALSE)


####################
# Pseudoreplicates #
####################

# consider everything with the same species id as a pseudoreplicate,
# combine E, and print out to a new file

# id <- c()
# id_process <- c()
# E_sch <- c()
# lnB0_sch <- c()
# E_D_sch <- c()
# T_h_sch <- c()
# T_pk_sch <- c()
# P_pk_sch <- c()
# habitat <- c()
# ConKingdom <- c()
# ConPhylum <- c()
# ConClass <- c()
# ConOrder <- c()
# ConFamily <- c()
# ConGenus <- c()
# ConSpecies <- c()
# TaxonomicGroup <- c()
# 
# sp_ids <- as.character(unique(autotroph_data$id_spp))
# 
# 
# for(i in 1:length(sp_ids)){
#   
#   data_subs <- autotroph_data[autotroph_data$id_spp == sp_ids[i],]
#   
#   id <- c(id, sp_ids[i])
#   id_process <- c(id_process, as.character(unique(data_subs$id_process)))
#   #E_sch <- c(E_sch, mean(data_subs$E_sch, na.rm = TRUE))
#   E_sch <- c(E_sch, weighted.mean(data_subs$E_sch, 1/(data_subs$E_SE_sch), na.rm = TRUE))
#   lnB0_sch <- c(lnB0_sch, mean(data_subs$lnB0_sch, na.rm = TRUE))
#   E_D_sch <- c(E_D_sch, mean(data_subs$E_D_sch, na.rm = TRUE))
#   T_h_sch <- c(T_h_sch, mean(data_subs$T_h_sch, na.rm = TRUE))
#   T_pk_sch <- c(T_pk_sch, mean(data_subs$T_pk_sch, na.rm = TRUE))
#   P_pk_sch <- c(P_pk_sch, mean(data_subs$P_pk_sch, na.rm = TRUE))
#   habitat <- c(habitat, unique(data_subs$habitat))
#   ConKingdom <- c(ConKingdom, unique(data_subs$ConKingdom))
#   ConPhylum <- c(ConPhylum, unique(data_subs$ConPhylum))
#   ConClass <- c(ConClass, unique(data_subs$ConClass))
#   ConOrder <- c(ConOrder, unique(data_subs$ConOrder))
#   ConFamily <- c(ConFamily, unique(data_subs$ConFamily))
#   ConGenus <- c(ConGenus, unique(data_subs$ConGenus))
#   ConSpecies <- c(ConSpecies, unique(data_subs$ConSpecies))
#   TaxonomicGroup <- c(TaxonomicGroup, unique(data_subs$TaxonomicGroup))
# }
# 
# autotroph_data_no_pseudo <- data.frame(id, id_process, E_sch, lnB0_sch, E_D_sch, T_h_sch, T_pk_sch, P_pk_sch,
#                                        habitat, ConKingdom, ConPhylum, ConClass, ConOrder, ConFamily, ConGenus, ConSpecies,
#                                        TaxonomicGroup)


# Alternative method: just take the best fit for each pseudoreplicate as representative
autotroph_data <- read.csv("../../Data/sofia_data/DataAutotrophs.csv")

autotroph_data_no_pseudo <- data.frame(matrix(NA, 0, ncol = ncol(autotroph_data)))
names(autotroph_data_no_pseudo) <- names(autotroph_data)

sp_ids <- as.character(unique(autotroph_data$id_spp))

for(i in 1:length(sp_ids)){
  
  subs_data <- autotroph_data[autotroph_data$id_spp == sp_ids[i],]
  
  if(length(subs_data$X) > 1){ # check for replicates
    
    # We're really only interested in E, so pick the fit with lowest error for E
    best_data <- subs_data[subs_data$E_std == min(subs_data$E_std, na.rm = TRUE),]
    
    autotroph_data_no_pseudo <- rbind(autotroph_data_no_pseudo, best_data)
    
  }
  else{ # if there were no replicates, just add it back to the dataset
    
    autotroph_data_no_pseudo <- rbind(autotroph_data_no_pseudo, subs_data)
    
  }
}

# add an extra column sort of for the python stuff later
autotroph_data_no_pseudo$Trophy <- "Autotroph"


# write to file
write.csv(autotroph_data_no_pseudo, "../../Data/sofia_data/DataAutotrophs_nopseudo.csv")

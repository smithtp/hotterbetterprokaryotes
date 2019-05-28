rm(list = ls())

require(ggplot2)
require(reshape)
require(grid)
require(gridExtra)

main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size=18, vjust=1),
        legend.text=element_text(size=14),
        strip.text.x = element_text(size = 16))

setwd("~/Documents/hotterbetterprokaryotes/Code/R_Scripts/")

data <- read.csv("../../Data/summaries/group_means_non_aggregated_taxonomic.csv")

# we need to make a new dataframe for each taxonomic group then bind them together I think

Kingdom_data <- data[data$ConKingdom_len > 4, ]
Kingdom_ES <- unique(Kingdom_data[, c('ConKingdom', 'ConKingdom', 'ConKingdom_ES', 'ConKingdom_ES_min', 'ConKingdom_ES_max', 'ConKingdom_len',
                                      "ConKingdom_ES_median")])
names(Kingdom_ES) <- c("Taxonomic_Group", "Kingdom", "E", "E_min", "E_max", "Length", "ES_median")
Kingdom_ES$Taxonomic_Group <- as.character(Kingdom_ES$Taxonomic_Group)
Kingdom_ES$Level <- "Kingdom"

Phylum_data <- data[data$ConPhylum_len > 4, ]
Phylum_data <- Phylum_data[!is.na(Phylum_data$ConPhylum),]
Phylum_ES <- unique(Phylum_data[, c('ConPhylum', 'ConKingdom', 'ConPhylum_ES', 'ConPhylum_ES_min', 'ConPhylum_ES_max', 'ConPhylum_len',
                                    "ConPhylum_ES_median")])
names(Phylum_ES) <- c("Taxonomic_Group", "Kingdom", "E", "E_min", "E_max", "Length", "ES_median")
Phylum_ES$Taxonomic_Group <- as.character(Phylum_ES$Taxonomic_Group)
Phylum_ES$Level <- "Phylum"

Class_data <- data[data$ConClass_len > 4, ]
Class_data <- Class_data[!is.na(Class_data$ConClass),]
Class_ES <- unique(Class_data[, c('ConClass', 'ConKingdom', 'ConClass_ES', 'ConClass_ES_min', 'ConClass_ES_max', 'ConClass_len',
                                  "ConClass_ES_median")])
names(Class_ES) <- c("Taxonomic_Group", "Kingdom", "E", "E_min", "E_max", "Length", "ES_median")
Class_ES$Taxonomic_Group <- as.character(Class_ES$Taxonomic_Group)
Class_ES$Level <- "Class"

Order_data <- data[data$ConOrder_len > 4, ]
Order_data <- Order_data[!is.na(Order_data$ConOrder),]
Order_ES <- unique(Order_data[, c('ConOrder', 'ConKingdom', 'ConOrder_ES', 'ConOrder_ES_min', 'ConOrder_ES_max', 'ConOrder_len', 
                                  "ConOrder_ES_median")])
names(Order_ES) <- c("Taxonomic_Group", "Kingdom", "E", "E_min", "E_max", "Length", "ES_median")
Order_ES$Taxonomic_Group <- as.character(Order_ES$Taxonomic_Group)
Order_ES$Level <- "Order"

combined_data <- rbind(Kingdom_ES, Phylum_ES, Class_ES, Order_ES)

# now a little addition for nicer y axis labels
combined_data$Taxonomic_Group_labels <- paste(combined_data$Taxonomic_Group, " (", combined_data$Length, ")", sep = "")

combined_data$Level <- factor(combined_data$Level, levels = c("Kingdom", "Phylum", "Class", "Order"))

main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 48),
        axis.text.y = element_text(size = 32),
        axis.title.y = element_text(size = 56),
        axis.title.x = element_text(size = 56),
        plot.title = element_text(size = 48, vjust=1),
        legend.text = element_text(size = 48),
        legend.title = element_text(size = 40),
        strip.text.x = element_text(size = 36),
        strip.text.y = element_text(size = 36))

png("../../Results/figures/ES_plot.png", height = 2400, width = 1600)
ES_plot <- ggplot(combined_data, aes(x = E, y = Taxonomic_Group_labels)) +
  geom_point(size = 10, aes(col = Kingdom)) +
  geom_errorbarh(aes(xmax=E_max, xmin=E_min, col = Kingdom), size = 3) +
  geom_point(aes(x = ES_median), size = 8, fill = "#CCCCCC", shape = 24) +
  geom_vline(aes(xintercept=0.65), linetype = 'dotted', size = 2) +
  main_theme +
  coord_cartesian(xlim=c(0, 3)) +
  ylab("Taxonomic Group") +
  xlab("Thermal Sensitivity (eV)") +
  facet_grid(Level ~ ., scales='free', space = "free")  +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(legend.title=element_blank(), legend.position=c(0.85,0.6)) +
  scale_color_manual(values=c("#FF9900", "#003399"))
ES_plot
dev.off()


############################################################
## Now repeat with functional groups rather than taxonomic #
############################################################

data_funct <- read.csv("../../Data/summaries/group_means_non_aggregated_funcgroups.csv")

# we need to make a new dataframe for each taxonomic group then bind them together I think

Pathogen_data <- data_funct[data_funct$Pathogen_len > 4 ,]
Pathogen_ES <- unique(Pathogen_data[, c('Pathogen', 'Pathogen_ES', 'Pathogen_ES_min', 'Pathogen_ES_max', 'Pathogen_len', 
                                        'Pathogen_ES_median')])
names(Pathogen_ES) <- c("Functional_Group", "E", "E_min", "E_max", "Length", "E_median")
Pathogen_ES$Functional_Group <- as.character(Pathogen_ES$Functional_Group)
Pathogen_ES$Level <- "Pathogen Status"

PathCarrier_data <- data_funct[data_funct$Pathogen_type_len > 4 ,]
PathCarrier_data <- PathCarrier_data[!is.na(PathCarrier_data$Pathogen_type),]
PathCarrier_ES <- unique(PathCarrier_data[, c('Pathogen_type', 'Pathogen_type_ES', 'Pathogen_type_ES_min', 'Pathogen_type_ES_max', 
                                              'Pathogen_type_len', 'Pathogen_type_ES_median')])
names(PathCarrier_ES) <- c("Functional_Group", "E", "E_min", "E_max", "Length", "E_median")
PathCarrier_ES$Functional_Group <- as.character(PathCarrier_ES$Functional_Group)
PathCarrier_ES$Level <- "Pathogen Carrier"
PathCarrier_ES <- PathCarrier_ES[PathCarrier_ES$Functional_Group != "Others",]

Metabolic_data <- data_funct[data_funct$Best_Guess_len > 4 ,]
Metabolic_ES <- unique(Metabolic_data[, c('Best_Guess', 'Best_Guess_ES', 'Best_Guess_ES_min', 'Best_Guess_ES_max', 'Best_Guess_len',
                                          'Best_Guess_ES_median')])
names(Metabolic_ES) <- c("Functional_Group", "E", "E_min", "E_max", "Length", "E_median")
Metabolic_ES$Functional_Group <- as.character(Metabolic_ES$Functional_Group)
# now we should put better names on these
Metabolic_ES$Functional_Group[1] <- "Acetogenesis"
Metabolic_ES$Functional_Group[2] <- "Aerobic Respiration"
Metabolic_ES$Functional_Group[3] <- "Fermentation"
Metabolic_ES$Functional_Group[4] <- "Photosynthesis"
Metabolic_ES$Functional_Group[5] <- "Sulfate Reduction"
Metabolic_ES$Functional_Group[6] <- "Methanogenesis"
Metabolic_ES$Functional_Group[7] <- "Sulfur Reduction"
Metabolic_ES$Functional_Group[8] <- "Hydrogen Oxidation"
Metabolic_ES$Functional_Group[9] <- "Iron Reduction"
Metabolic_ES$Functional_Group[10] <- "Sulfur Oxidation"

Metabolic_ES$Level <- "Energy Generation"


TempPref_data <- data_funct[data_funct$TempPref_len > 4 ,]
TempPref_ES <- unique(TempPref_data[, c('TempPref', 'TempPref_ES', 'TempPref_ES_min', 'TempPref_ES_max', 'TempPref_len', 
                                        'TempPref_ES_median')])
TempPref_ES <- TempPref_ES[!is.na(TempPref_ES$TempPref),]
names(TempPref_ES) <- c("Functional_Group", "E", "E_min", "E_max", "Length", "E_median")
TempPref_ES$Functional_Group <- as.character(TempPref_ES$Functional_Group)
TempPref_ES$Level <- "Temperature Niche"

combined_data_funct <- rbind(Pathogen_ES, PathCarrier_ES, Metabolic_ES, TempPref_ES)

# now a little addition for nicer y axis labels
combined_data_funct$Functional_Group_labels <- paste(combined_data_funct$Functional_Group, " (", combined_data_funct$Length, ")", sep = "")

combined_data_funct$Level <- factor(combined_data_funct$Level, levels = c("Pathogen Status", "Pathogen Carrier", "Energy Generation", "Temperature Niche"))

png("../../Results/figures/ES_plot_funct.png", height = 1000, width = 1600)
ES_plot_funct <- ggplot(combined_data_funct, aes(x = E, y = Functional_Group_labels)) +
  geom_point(size = 10) +
  geom_errorbarh(aes(xmax=E_max, xmin=E_min), size = 3) +
  geom_point(aes(x = E_median), size = 8, fill = "#CCCCCC", shape= 24) +
  geom_vline(aes(xintercept=0.65), linetype = 'dotted', size = 2) +
  main_theme +
  coord_cartesian(xlim=c(0, 3)) +
  ylab("Functional Group") +
  xlab("Thermal Sensitivity (eV)") +
  facet_grid(Level ~ ., scales='free', space = "free")  +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(legend.title=element_blank(), legend.position=c(0.85,0.6))
ES_plot_funct
dev.off()



############################################################
## Now repeat for the autotroph dataset                    #
############################################################

data_autotrophs <- read.csv("../../Data/sofia_data/DataAutotrophs_nopseudo.csv")

# we're going to need to do the bootstrap of weighted mean on E for different groups

# implement bootstrapping code for CIs, like the python code

resample_all <- function(x){
  boot_nums = c()
  len_x = length(x)
  for (i in 1:1000){
    boot_nums <- c(boot_nums, mean(sample(x, len_x, replace=TRUE)))
  }
  return( boot_nums )
} 

bootstrap_upper <- function(x){
  sample <- resample_all(x)
  quant = as.single(quantile(sample, 0.975))[1]
  return(quant)
}

bootstrap_mean <- function(x){
  sample <- resample_all(x)
  mn <- as.single(mean(sample))[1]
  return(mn)
}

bootstrap_lower <- function(x){
  sample <- resample_all(x)
  quant = as.single(quantile(sample, 0.025))[1]
  return(quant)
}

# make a dataframe to put things in
combined_data_autotroph <- data.frame(matrix(NA, 7, 8))
names(combined_data_autotroph) <- c("Group", "E", "E_min", "E_max", "Length", "E_median", "Level", "Group_labels")

# now bootstrap E first for terrestrial/aqauatic
terr_data <- data_autotrophs[data_autotrophs$habitat == "Terrestrial",]
combined_data_autotroph$Group[1] <- "Terrestrial"
combined_data_autotroph$E[1] <- bootstrap_mean(terr_data$E)
combined_data_autotroph$E_min[1] <- bootstrap_lower(terr_data$E)
combined_data_autotroph$E_max[1] <- bootstrap_upper(terr_data$E)
combined_data_autotroph$Length[1] <- length(terr_data$E)
combined_data_autotroph$E_median[1] <- median(terr_data$E)
combined_data_autotroph$Level[1] <- "Habitat"

aq_data <- data_autotrophs[data_autotrophs$habitat == "Aquatic",]
combined_data_autotroph$Group[2] <- "Aquatic"
combined_data_autotroph$E[2] <- bootstrap_mean(aq_data$E)
combined_data_autotroph$E_min[2] <- bootstrap_lower(aq_data$E)
combined_data_autotroph$E_max[2] <- bootstrap_upper(aq_data$E)
combined_data_autotroph$Length[2] <- length(aq_data$E)
combined_data_autotroph$E_median[2] <- median(aq_data$E)
combined_data_autotroph$Level[2] <- "Habitat"

# now do the taxonomic group levels
plants <- data_autotrophs[data_autotrophs$TaxonomicGroup == "Vascular Plants",]
combined_data_autotroph$Group[3] <- "Vascular Plants"
combined_data_autotroph$E[3] <- bootstrap_mean(plants$E)
combined_data_autotroph$E_min[3] <- bootstrap_lower(plants$E)
combined_data_autotroph$E_max[3] <- bootstrap_upper(plants$E)
combined_data_autotroph$Length[3] <- length(plants$E)
combined_data_autotroph$E_median[3] <- median(plants$E)
combined_data_autotroph$Level[3] <- "Taxonomic Group"

greenalgae <- data_autotrophs[data_autotrophs$TaxonomicGroup == "Green Algae",]
combined_data_autotroph$Group[4] <- "Green Algae"
combined_data_autotroph$E[4] <- bootstrap_mean(greenalgae$E)
combined_data_autotroph$E_min[4] <- bootstrap_lower(greenalgae$E)
combined_data_autotroph$E_max[4] <- bootstrap_upper(greenalgae$E)
combined_data_autotroph$Length[4] <- length(greenalgae$E)
combined_data_autotroph$E_median[4] <- median(greenalgae$E)
combined_data_autotroph$Level[4] <- "Taxonomic Group"

redalgae <- data_autotrophs[data_autotrophs$TaxonomicGroup == "Red Algae",]
combined_data_autotroph$Group[5] <- "Red Algae"
combined_data_autotroph$E[5] <- bootstrap_mean(redalgae$E)
combined_data_autotroph$E_min[5] <- bootstrap_lower(redalgae$E)
combined_data_autotroph$E_max[5] <- bootstrap_upper(redalgae$E)
combined_data_autotroph$Length[5] <- length(redalgae$E)
combined_data_autotroph$E_median[5] <- median(redalgae$E)
combined_data_autotroph$Level[5] <- "Taxonomic Group"

brownalgae <- data_autotrophs[data_autotrophs$TaxonomicGroup == "Brown Algae",]
combined_data_autotroph$Group[6] <- "Brown Algae"
combined_data_autotroph$E[6] <- bootstrap_mean(brownalgae$E)
combined_data_autotroph$E_min[6] <- bootstrap_lower(brownalgae$E)
combined_data_autotroph$E_max[6] <- bootstrap_upper(brownalgae$E)
combined_data_autotroph$Length[6] <- length(brownalgae$E)
combined_data_autotroph$E_median[6] <- median(brownalgae$E)
combined_data_autotroph$Level[6] <- "Taxonomic Group"

mosses <- data_autotrophs[data_autotrophs$TaxonomicGroup == "Mosses",]
combined_data_autotroph$Group[7] <- "Mosses"
combined_data_autotroph$E[7] <- bootstrap_mean(mosses$E)
combined_data_autotroph$E_min[7] <- bootstrap_lower(mosses$E)
combined_data_autotroph$E_max[7] <- bootstrap_upper(mosses$E)
combined_data_autotroph$Length[7] <- length(mosses$E)
combined_data_autotroph$E_median[7] <- median(mosses$E)
combined_data_autotroph$Level[7] <- "Taxonomic Group"

combined_data_autotroph$Group_labels <- paste(combined_data_autotroph$Group, " (", combined_data_autotroph$Length, ")", sep = "")

png("../../Results/figures/ES_plot_autotroph.png", height = 600, width = 1600)
ES_plot_autotroph <- ggplot(combined_data_autotroph, aes(x = E, y = Group_labels)) +
  geom_point(size = 10) +
  geom_errorbarh(aes(xmax=E_max, xmin=E_min), size = 3) +
  geom_point(aes(x = E_median), size = 8, fill = "#CCCCCC", shape= 24) +
  geom_vline(aes(xintercept=0.65), linetype = 'dotted', size = 2) +
  main_theme +
  #coord_cartesian(xlim=c(0, 3)) +
  #ylab("Group") +
  xlab("Thermal Sensitivity (eV)") +
  facet_grid(Level ~ ., scales='free', space = "free")  +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(legend.title=element_blank(), legend.position=c(0.85,0.6),
        axis.title.y = element_blank())
ES_plot_autotroph
dev.off()

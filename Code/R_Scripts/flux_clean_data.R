# tidy up flux data

rm(list = ls())

setwd("~/Documents/hotterbetterprokaryotes/Data/summaries/")

raw_data <- read.csv("non_aggregate_data_fluxes.csv")

best_data <- raw_data[raw_data$Trait %in% c(
  "Caffeine Degradaded Per second (g)",
  "CH4 produced",
  "Dehalogenation",
  "Fe2+ Oxidation Rate",
  "Fe2+ Oxidised",
  "Fe3+ Oxidation Rate",
  "H2S Production",
  "Methanogenesis",
  "mM Fe(II) production",
  "Nitrate Production (mM)",
  "Nitrogen Removal Rate",
  "SO4 Reduction Rate",
  "Sulfate Reduction Rate",
  "Sulphide Production",
  "Sulphur Oxidation Rate",
  "Time To Double Fe3+ concentration",
  "Time To Increase SO4 to 3000mg/L",
  "Time To oxidise 2000mg/L Fe2+",
  "uMols methane"),]

write.csv(best_data, "non_aggregate_data_fluxes_best.csv", row.names = FALSE)

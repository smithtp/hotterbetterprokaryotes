## Re-attach the temperature preferences based on segmented analysis

rm(list = ls())

setwd("~/Documents/hotterbetterprokaryotes/Code/R_Scripts/")

data <- read.csv("../../Data/summaries/non_aggregate_data.csv")

names(data)

# ggplot(data, aes(x = Est_Tpk, y = Est_Response, col = TempPref)) +
#   geom_point()

#write a new ifelse statement here

data$TempPref <- NA

for (i in 1:length(data$Est_Tpk)){
  if (!is.na(data$Est_Tpk[i])){
    if (data$ConKingdom[i] == "Bacteria"){
      if (data$Est_Tpk[i] > 313.65){ #40.5C
        data$TempPref[i] <- "Thermophile"
      }
      else{
        data$TempPref[i] <- "Mesophile"
      }
    }
    else if (data$ConKingdom[i] == "Archaea"){
      if (data$Est_Tpk[i] > 319.15){ # 46C
        data$TempPref[i] <- "Thermophile"
      }
      else{
        data$TempPref[i] <- "Mesophile"
      }
    }
    else{
      data$TempPref[i] <- "NA"
    }
  }
  else{
    data$TempPref[i] <- "NA"
  }
}

# ggplot(data, aes(x = Est_Tpk, y = Est_Response, col = TempPref)) +
#   geom_point()

write.csv(data, "../../Data/summaries/non_aggregate_data_fix.csv", row.names = FALSE)

bacteria <- data[data$ConKingdom == "Bacteria",]

write.csv(bacteria, "../../Data/summaries/non_aggregate_data_fix_bacteria.csv", row.names = FALSE)

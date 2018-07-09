## Using "segmented" R package to find analysis break-points in the data

rm(list = ls())

#install.packages("segmented")
library("segmented")
library("ggplot2")

setwd("~/Documents/hotterbetterprokaryotes/Code/R_Scripts/")

data <- read.csv("../../Data/summaries/non_aggregate_data.csv")

bacteria_data <- data[data$ConKingdom == "Bacteria",]
archaea_data <- data[data$ConKingdom == "Archaea",]

ggplot(data, aes(x = 1/Est_Tpk, y = log(Est_Response), col = ConKingdom)) +
  geom_point()

ggplot(bacteria_data, aes(x = 1/Est_Tpk, y = log(Est_Response))) +
  geom_point()

ggplot(archaea_data, aes(x = 1/Est_Tpk, y = log(Est_Response))) +
  geom_point()

bacteria <- bacteria_data

# lets just use a linear form of BA model

k <- 8.62e-5 

x = 1 / (k*bacteria$Est_Tpk)
y = log(bacteria$Est_Response)

boltzmann_linear <- lm(y ~ x)

os <- segmented(boltzmann_linear, seg.Z =~ x) # break point at 36.97?!
plot(x,y)
plot(os, add=T)
summary.segmented(os)

os$psi[2] # <-- this is the break point

bacteria_breakpoint_celsius <- (1/(k*os$psi[2]))-273.15
bacteria_breakpoint_celsius

# repeat for archaea

x = 1 / (k*archaea_data$Est_Tpk)
y = log(archaea_data$Est_Response)

boltzmann_linear <- lm(y ~ x)

os <- segmented(boltzmann_linear, seg.Z =~ x) # break point at 36.97?!
plot.segmented(os)
summary.segmented(os)

os$psi[2] # <-- this is the break point

archaea_breakpoint_celsius <- (1/(k*os$psi[2]))-273.15
archaea_breakpoint_celsius

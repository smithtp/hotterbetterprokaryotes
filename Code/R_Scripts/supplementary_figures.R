## supplementary data figures

rm(list = ls())

library("ggplot2")

setwd("~/Documents/hotterbetterprokaryotes/Code/R_Scripts/")

main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        plot.title = element_text(size=18, vjust=1),
        legend.text=element_text(size=18),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 16))

original_data <- read.csv("../../Data/database.csv")
original_growth_data <- original_data[original_data$StandardisedTraitName == "Specific Growth Rate",]

growth_rates_plot <- ggplot(original_growth_data, aes(x = ConTemp, y = StandardisedTraitValue*3600, col = ConKingdom)) +
    geom_point(alpha = 0.8) +
    ylab("Specific Growth Rate (per hour)") +
    xlab("Growth Temperature (Celsius)") +
    main_theme +
    theme(legend.title=element_blank(), legend.position=c(0.1,0.9)) +
    guides(colour = guide_legend(override.aes = list(size=5))) # increase legend point size
growth_rates_plot

ggsave(file = '../../Results/supplement_figs/all_growth_rates.png', growth_rates_plot, width = 12, height = 8)

fitting_data <- read.csv("../../Data/summaries/summary.csv")

### -- show Tpk with lab growth temp -- ### (microbes grow best near the temp they're routinely grown at!)
### though whether this is due to fast adaptation to lab conditions or species sorting during isolation is a point for discussion

Tpk_vs_labtemp <- ggplot(fitting_data, aes(x = ConLabGrowthTemp, y = Est_Tpk-273.15, col = ConKingdom)) +
      geom_point(size = 2, alpha = 0.8) +
      geom_abline(aes(slope = 1, intercept = 0)) +
      #geom_smooth(method =lm, se = FALSE) +
      xlim(0, 110) +
      ylim(0, 110) +
      geom_smooth(method = lm) +
      main_theme +
      theme(legend.title=element_blank(), legend.position=c(0.1,0.9)) +
      xlab("Laboratory Growth Temperature (Celsius)") +
      ylab("Peak Growth Temperature (Celsius)") +
      guides(colour = guide_legend(override.aes = list(size=5)))
Tpk_vs_labtemp

ggsave(file = '../../Results/supplement_figs/lab_temp_vs_peak_temp.png', Tpk_vs_labtemp, width = 12, height = 8)

# actually do the linear models
fitting_data$Tpk_C <- fitting_data$Est_Tpk-273.15

bact_lm <- lm(Tpk_C ~ ConLabGrowthTemp, data = fitting_data[fitting_data$ConKingdom == "Bacteria",])
summary(bact_lm)
confint(bact_lm, "ConLabGrowthTemp", level=0.95)

arch_lm <- lm(Tpk_C ~ ConLabGrowthTemp, data = fitting_data[fitting_data$ConKingdom == "Archaea",])
summary(arch_lm)
confint(arch_lm, "ConLabGrowthTemp", level=0.95)

### -- show that E doesn't rise with Tpk -- ### (could confound analysis) if increase in specialists with temp
### but what about E with max rate?
# E vs Tpk

E_vs_Tpk <- ggplot(fitting_data, aes(x = Est_Tpk-273.15, y = E)) +
            geom_point() +
            ylim(0, 5) +
            xlab("Peak Growth Temperature") +
            ylab("Activation Energy")

ggsave(file = '../../Results/supplement_figs/E_vs_Tpk.png', E_vs_Tpk, width = 12, height = 8)

# now also test the corellation formally:

cor.test(fitting_data$Est_Tpk, fitting_data$E)

growth_data <- fitting_data[fitting_data$Trait == "Specific Growth Rate",]

ggplot(growth_data, aes(x = E, y = Max_response)) +
  geom_point() +
  xlim(0, 5)

### -- look at parameter errors with temperature -- ## (if errors rise with temp, can see why weighted BA drops E compared to normal?)
# E error vs Tpk

ggplot(fitting_data, aes(x = Est_Tpk, y = E_std)) +
  geom_point() +
  ylim(0, 5)

# Tpk error vs Tpk

ggplot(fitting_data, aes(x = Est_Tpk, y = Tpk_std)) +
  geom_point() +
  ylim(0, 10)

# Ppk error vs Tpk

ggplot(growth_data, aes(x = Est_Tpk, y = Response_std)) +
  geom_point() +
  ylim(0, 20) 
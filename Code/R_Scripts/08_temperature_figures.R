rm(list = ls())
library("ggplot2")
library("minpack.lm")
library("gtable")
library("grid")
library("gridExtra")

setwd("~/Documents/hotterbetterprokaryotes/Code/R_Scripts/")

main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size=20, vjust=1),
        legend.text=element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

data <- read.csv("../../Data/summaries/non_aggregate_data_fix.csv")

bacteria <- data[data$ConKingdom == "Bacteria",]
bacteria <- bacteria[!is.na(bacteria$TempPref),]
archaea <- data[data$ConKingdom == "Archaea",]
archaea <- archaea[!is.na(archaea$TempPref),]

arrh_data <- read.csv("../../Data/summaries/all_means_non_aggregated.csv")

###################### Boltzmann - Arrhenius model.
Boltzmann.Arrhenius <- function(lnB0, E, temp) {
  
  # Boltzmann's constant. Units imply that E is in eV.
  k <- 8.62e-5 
  
  # lnB0 is the normalization constant.  
  # E is the activation energy.
  # Tref is the standardization temperature (in K).
  
  calc <- lnB0 - E/k * (1/temp)
  
  return(calc)
}

###### Super fancy plotting function

fancy_plot <- function(p1, p2){
  ## extract gtable
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))
  
  ## overlap the panel of the 2nd plot on that of the 1st plot
  pp <- c(subset(g1$layout, name=="panel", se=t:r))
  
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b, 
                       pp$l)
  
  g <- gtable_add_grob(g1, g1$grobs[[which(g1$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)
  #g <- gtable_add_grob(g1, g1$grobs[[which(g1$layout$name=="guide-box")]], pp$t, pp$l, pp$b, pp$l) # this line for the legend
  
  ## steal axis from second plot and modify
  ia <- which(g2$layout$name == "axis-b")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[[2]]
  
  ## switch position of ticks and labels
  ax$heights <- rev(ax$heights)
  ax$grobs <- rev(ax$grobs)
  #ax$grobs[[2]]$y <- ax$grobs[[2]]$y - unit(1, "npc") + unit(0.15, "cm") # remove this line when removing axis ticks!
  
  ## modify existing row to be tall enough for axis
  g$heights[[2]] <- g$heights[g2$layout[ia,]$t]
  
  ## add new axis
  g <- gtable_add_grob(g, ax, 2, 4, 2, 4)
  
  ## add new row for upper axis label
  #g <- gtable_add_rows(g, g2$heights[1], 1) ## removed this because its not placing properly, will add label manually later
  #g <- gtable_add_grob(g, g2$grob[[12]], 2, 4, 2, 4)
  
  # draw it
  grid.draw(g)
}

### new and improved with legend carry-over!

fancy_plot_legend <- function(p1, p2){
  ## extract gtable
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))
  
  ## overlap the panel of the 2nd plot on that of the 1st plot
  pp <- c(subset(g1$layout, name=="panel", se=t:r))
  
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b, 
                       pp$l)
  
  g <- gtable_add_grob(g1, g1$grobs[[which(g1$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)
  g <- gtable_add_grob(g1, g1$grobs[[which(g1$layout$name=="guide-box")]], pp$t, pp$l, pp$b, pp$l) # this line for the legend
  
  ## steal axis from second plot and modify
  ia <- which(g2$layout$name == "axis-b")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[[2]]
  
  ## switch position of ticks and labels
  ax$heights <- rev(ax$heights)
  ax$grobs <- rev(ax$grobs)
  #ax$grobs[[2]]$y <- ax$grobs[[2]]$y - unit(1, "npc") + unit(0.15, "cm") # remove this line when removing axis ticks!
  
  ## modify existing row to be tall enough for axis
  g$heights[[2]] <- g$heights[g2$layout[ia,]$t]
  
  ## add new axis
  g <- gtable_add_grob(g, ax, 2, 4, 2, 4)
  
  ## add new row for upper axis label
  #g <- gtable_add_rows(g, g2$heights[1], 1) # 5.5pt btw, hence changed in below line as need more space
  g <- gtable_add_rows(g, unit(20, "pt"), 1)
  g <- gtable_add_grob(g, g2$grob[[12]], 2, 4, 2, 4)
  
  # draw it
  grid.draw(g)
}

########
# generate the bacteria BA curve

Kingdom_BA <- unique(arrh_data[, c('ConKingdom', 'ConKingdom', 'ConKingdom_EG', 'ConKingdom_EG_min', 'ConKingdom_EG_max', 'ConKingdom_B0')])
bacteria_BA <- Kingdom_BA[Kingdom_BA$ConKingdom == "Bacteria",]
archaea_BA <- Kingdom_BA[Kingdom_BA$ConKingdom == "Archaea",]

# generate a curve from the resultant model:
bacteria_temps <- seq(273.15, 373.15, 0.1)
bacteria_responses <- Boltzmann.Arrhenius(bacteria_BA$ConKingdom_B0, bacteria_BA$ConKingdom_EG, bacteria_temps)
bacteria_curve_data <- data.frame(bacteria_temps, bacteria_responses)

k=8.62e-5
ggplot(bacteria, aes(x = 1/(k*Est_Tpk), y = log(Est_Response))) +
  geom_point() +
  geom_line(data = bacteria_curve_data, aes(x = 1/(k*bacteria_temps), y = bacteria_responses))

# with 2 X-axes?
# first a quick look at how celsius and 1/kT scale
celsius <- seq(0, 100, 0.5)
kT <- 1/(k*(celsius+273.15))
temperatureframe <- data.frame(celsius, kT)

# now figure out the plotting
p1 <- ggplot(bacteria, aes(x = 1/(k*Est_Tpk), y = log(Est_Response))) +
  geom_point() +
  geom_line(data = bacteria_curve_data, aes(x = 1/(k*bacteria_temps), y = bacteria_responses)) +
  scale_x_continuous(name="1/ kT",limits=c(31,43),
                     breaks=c(32.5, 35, 37.5, 40, 42.5))

p2 <- ggplot(bacteria, aes(x = 1/(k*Est_Tpk), y = log(Est_Response))) +
  geom_point() +
  geom_line(data = bacteria_curve_data, aes(x = 1/(k*bacteria_temps), y = bacteria_responses)) +
  scale_x_continuous(name="Temperature (celsius)",limits=c(31,42),
                     breaks=c(42.47,39.57,37.05,34.82,32.85,31.09),
                     labels=c(0,20,40,60,80,100) ## labels convert to celsius
  )

fancy_plot(p1, p2)

## look at archaea too...

archaea_temps <- seq(273.15, 375.15, 0.1)
archaea_responses <- Boltzmann.Arrhenius(archaea_BA$ConKingdom_B0, archaea_BA$ConKingdom_EG, archaea_temps)
archaea_curve_data <- data.frame(archaea_temps, archaea_responses)

k=8.62e-5
ggplot(archaea, aes(x = 1/(k*Est_Tpk), y = log(Est_Response))) +
  geom_point() +
  geom_line(data = archaea_curve_data, aes(x = 1/(k*archaea_temps), y = archaea_responses))


p1 <- ggplot(archaea, aes(x = 1/(k*Est_Tpk), y = log(Est_Response))) +
  geom_point() +
  geom_line(data = archaea_curve_data, aes(x = 1/(k*archaea_temps), y = archaea_responses)) +
  scale_x_continuous(name="1/ kT",limits=c(31,43),
                     breaks=c(32.5, 35, 37.5, 40, 42.5))

p2 <- ggplot(archaea, aes(x = 1/(k*Est_Tpk), y = log(Est_Response))) +
  geom_point() +
  geom_line(data = archaea_curve_data, aes(x = 1/(k*archaea_temps), y = archaea_responses)) +
  scale_x_continuous(name="Temperature (celsius)",limits=c(31,42),
                     breaks=c(42.47,39.57,37.05,34.82,32.85,31.09),
                     labels=c(0,20,40,60,80,100)) ## labels convert to celsius

fancy_plot(p1, p2)                     
      

###################################################################
################## MESOPHILE vs THERMOPHILE #######################
###################################################################

TempPref_BA <- unique(arrh_data[, c('ConKingdom', 'TempPref', 'TempPref_EG', 'TempPref_EG_min', 'TempPref_EG_max', 'TempPref_B0')])
TempPref_BA <- TempPref_BA[!is.na(TempPref_BA$TempPref),] # remove weird NAs

### bacteria first

bact_meso <- bacteria[bacteria$TempPref == "Mesophile",]
bact_thermo <- bacteria[bacteria$TempPref == "Thermophile",]

TempPref_BA_bacteria <- TempPref_BA[TempPref_BA$ConKingdom == "Bacteria",]
bacteria_mesoBA <- TempPref_BA_bacteria[TempPref_BA_bacteria$TempPref == "Mesophile",]
bacteria_thermoBA <- TempPref_BA_bacteria[TempPref_BA_bacteria$TempPref == "Thermophile",]

# should generate some lines for mesophile and thermophile models first
meso_temps <- seq(273.15, 313.65, 0.1)
meso_responses <- Boltzmann.Arrhenius(bacteria_mesoBA$TempPref_B0, bacteria_mesoBA$TempPref_EG, meso_temps)
meso_curve <- data.frame(meso_temps, meso_responses)
meso_curve$TempPref <- "Mesophile"

thermo_temps <- seq(313.65, 362, 0.1)
thermo_responses <- Boltzmann.Arrhenius(bacteria_thermoBA$TempPref_B0, bacteria_thermoBA$TempPref_EG, thermo_temps)
thermo_curve <- data.frame(thermo_temps, thermo_responses)
thermo_curve$TempPref <- "Thermophile"

p1 <- ggplot(bacteria, aes(x = 1/(k*Est_Tpk), y = log(Est_Response), col = TempPref)) +
  geom_point(size =3, alpha = 0.8) +
  geom_line(data = meso_curve, aes(x = 1/(k*meso_temps), y = meso_responses), size = 1.5) +
  geom_line(data = thermo_curve, aes(x = 1/(k*thermo_temps), y = thermo_responses), size = 1.5) +
  #xlab(expression("1/"~italic("kT"))) +
  ylab("log(maximum growth rate, per day)") +
  scale_x_continuous(name=expression("1/"~italic("kT"["pk"])),limits=c(31.5,41.5),
                   breaks=c(32.5, 35, 37.5, 40)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.justification=c(1.1,1.1), legend.position=c(1,1))

p2 <- ggplot(meso_curve, aes(x = 1/(k*meso_temps), y = meso_responses)) +
  geom_line() +
  xlab("Temperature (celsius)") +
  scale_x_continuous(limits=c(31.5,41.5),
                     breaks=c(39.57,37.05,34.82,32.85),
                     labels=c(20,40,60,80)) + ## labels convert to celsius
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 18, vjust =1))

  
png("../../Results/figures/bacteria_TempPref_weighted.png", width = 720, height = 480)
fancy_plot_legend(p1, p2)
dev.off()

## and an extra plot on a celsius axis for posters/presentations?

png("../../Results/figures/bacteria_TempPref_weighted_celsius.png", width = 720, height = 480)
ggplot(bacteria, aes(x = Est_Tpk-273.15, y = log(Est_Response), col = TempPref)) +
  geom_point(size =3, alpha = 0.8) +
  geom_line(data = meso_curve, aes(x = meso_temps-273.15, y = meso_responses), size = 1.5) +
  geom_line(data = thermo_curve, aes(x = thermo_temps-273.15, y = thermo_responses), size = 1.5) +
  ylab("log(maximum growth rate, per day)") +
  scale_x_continuous(name="Temperature (Celsius)",limits=c(0, 90)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.justification=c(1.1,1.1), legend.position=c(1,0.3))
dev.off()

## another alternative could be the arrhenius plot with flipped axes:

p1 <- ggplot(bacteria, aes(x = 1/(k*Est_Tpk), y = log(Est_Response), col = TempPref)) +
  geom_point(size =3, alpha = 0.8) +
  geom_line(data = meso_curve, aes(x = 1/(k*meso_temps), y = meso_responses), size = 1.5) +
  geom_line(data = thermo_curve, aes(x = 1/(k*thermo_temps), y = thermo_responses), size = 1.5) +
  ylab("log(maximum growth rate, per day)") +
  scale_x_continuous(trans = "reverse", name=expression("1/"~italic("kT"["pk"])),limits=c(41.5,31.5),
                     breaks=c(32.5, 35, 37.5, 40)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.justification=c(1.1,1.1), legend.position=c(1,0.3))

p2 <- ggplot(meso_curve, aes(x = 1/(k*meso_temps), y = meso_responses)) +
  geom_line() +
  xlab("Temperature (celsius)") +
  scale_x_continuous(trans = "reverse", limits=c(41.5,31.5),
                     breaks=c(39.57,37.05,34.82,32.85),
                     labels=c(20,40,60,80)) + ## labels convert to celsius
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 18, vjust =1))


png("../../Results/figures/bacteria_TempPref_weighted_axis_flip.png", width = 720, height = 480)
fancy_plot_legend(p1, p2)
dev.off()

############ now repeat for archaea... ###########

arch_meso <- archaea[archaea$TempPref == "Mesophile",]
arch_thermo <- archaea[archaea$TempPref == "Thermophile",]

TempPref_BA_archaea <- TempPref_BA[TempPref_BA$ConKingdom == "Archaea",]
archaea_mesoBA <- TempPref_BA_archaea[TempPref_BA_archaea$TempPref == "Mesophile",]
archaea_thermoBA <- TempPref_BA_archaea[TempPref_BA_archaea$TempPref == "Thermophile",]

# should generate some lines for mesophile and thermophile models first
meso_temps <- seq(286, 319.15, 0.1)
meso_responses <- Boltzmann.Arrhenius(archaea_mesoBA$TempPref_B0, archaea_mesoBA$TempPref_EG, meso_temps)
meso_curve <- data.frame(meso_temps, meso_responses)
meso_curve$TempPref <- "Mesophile"

thermo_temps <- seq(322.79, 377.3, 0.1)
thermo_responses <- Boltzmann.Arrhenius(archaea_thermoBA$TempPref_B0, archaea_thermoBA$TempPref_EG, thermo_temps)
thermo_curve <- data.frame(thermo_temps, thermo_responses)
thermo_curve$TempPref <- "Thermophile"

p1 <- ggplot(archaea, aes(x = 1/(k*Est_Tpk), y = log(Est_Response), col = TempPref)) +
  geom_point(size = 3, alpha = 0.8) +
  #geom_errorbar(aes(x = 1/(k*Est_Tpk), ymax = log(Est_Response + Est_Response_std), ymin = log(Est_Response - Est_Response_std))) +
  geom_line(data = meso_curve, aes(x = 1/(k*meso_temps), y = meso_responses), size = 1.5) +
  geom_line(data = thermo_curve, aes(x = 1/(k*thermo_temps), y = thermo_responses), size = 1.5) +
  ylab("log(maximum growth rate, per day)") +
  scale_y_continuous(limits=c(-2,4)) +
  scale_x_continuous(name=expression("1/"~italic("kT"["pk"])),limits=c(30,41),
                     breaks=c(30, 32.5, 35, 37.5, 40)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.justification=c(1.1,1.1), legend.position=c(1,1))

p2 <- ggplot(meso_curve, aes(x = 1/(k*meso_temps), y = meso_responses)) +
  geom_line() +
  scale_x_continuous(name="Temperature (celsius)",limits=c(30,41),
                     breaks=c(39.57,37.05,34.82,32.85,31.09),
                     labels=c(20,40,60,80,100)) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 18, vjust =1))

png("../../Results/figures/archaea_TempPref_weighted.png", width = 720, height = 480)
fancy_plot_legend(p1, p2)
dev.off()

# celsius

png("../../Results/figures/archaea_TempPref_weighted_celsius.png", width = 720, height = 480)
ggplot(archaea, aes(x = Est_Tpk-273.15, y = log(Est_Response), col = TempPref)) +
  geom_point(size =3, alpha = 0.8) +
  geom_line(data = meso_curve, aes(x = meso_temps-273.15, y = meso_responses), size = 1.5) +
  geom_line(data = thermo_curve, aes(x = thermo_temps-273.15, y = thermo_responses), size = 1.5) +
  xlab(expression("1/"~italic("kT"))) +
  ylab("log(maximum growth rate, per day)") +
  scale_x_continuous(name="Temperature (Celsius)",limits=c(0, 110)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.justification=c(1.1,1.1), legend.position=c(0.3,1))
dev.off()

## flipped axes:

p1 <- ggplot(archaea, aes(x = 1/(k*Est_Tpk), y = log(Est_Response), col = TempPref)) +
  geom_point(size =3, alpha = 0.8) +
  geom_line(data = meso_curve, aes(x = 1/(k*meso_temps), y = meso_responses), size = 1.5) +
  geom_line(data = thermo_curve, aes(x = 1/(k*thermo_temps), y = thermo_responses), size = 1.5) +
  ylab("log(maximum growth rate, per day)") +
  scale_x_continuous(trans = "reverse", name=expression("1/"~italic("kT"["pk"])),limits=c(41,30),
                     breaks=c(30, 32.5, 35, 37.5, 40)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D")) +
  scale_y_continuous(limits=c(-2,4)) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.justification=c(1.1,1.1), legend.position=c(0.3,1))

p2 <- ggplot(meso_curve, aes(x = 1/(k*meso_temps), y = meso_responses)) +
  geom_line() +
  xlab("Temperature (celsius)") +
  scale_x_continuous(trans = "reverse", limits=c(41,30),
                     breaks=c(39.57,37.05,34.82,32.85,31.09),
                     labels=c(20,40,60,80,100)) + ## labels convert to celsius
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 18, vjust =1))


png("../../Results/figures/archaea_TempPref_weighted_axis_flip.png", width = 720, height = 480)
fancy_plot_legend(p1, p2)
dev.off()

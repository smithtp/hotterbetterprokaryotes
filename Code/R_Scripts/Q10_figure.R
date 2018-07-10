####################################################################
# Use warming example calculation to produce a figure              #
# for all autotroph/heterotroph, prokaryote/eukaryote combinations #
####################################################################
library("ggplot2")
library("directlabels")

rm(list = ls())
graphics.off()

setwd("~/Documents/hotterbetterprokaryotes/Code/R_Scripts/")

main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        plot.title = element_text(size=18, vjust=1),
        legend.text=element_text(size=18),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16))

# set parameter values
k <- 8.617e-5 # boltzmann constant
Ef <- 0.65 # presumed E for fungi and other eukaryotic heterotrophs
Ea <- 0.65 # autotroph E (see fig 3)
Ep <- 0.87 # average intraspecific E for mesophilic bacteria
T1 <- 293.15 # 20C
T2 <- T1+10 # 10C change
C1 <- 2/(exp(-Ef/(k*T1))) # constant for eukaryotes/fungi
C2 <- 2/(exp(-Ep/(k*T1))) # constant for prokaryotes
C3 <- 2/(exp(-Ea/(k*T1))) # constant for autotrophs

# standard Q10 with 0.65eV for comparison
Q10_standard <- (C1*exp(-Ef/(k*T2)))/(C1*exp(-Ef/(k*T1)))


Q10 <- c()
change <- c()
auto_hetero <- c()
prok_euk <- c()
eco_E <- c()


for(i in 1:100){
  
  for(j in 1:100){
    
    d <- 0.01*i # auto vs hetero
    y <- 0.01*j # prok vs euk
    
    auto_hetero <- c(auto_hetero, d)
    prok_euk <- c(prok_euk, y)
    
    # first part autotroph vs heterotroph, second part bacteria vs eukaryotes within heterotrophs
    N1 <- d*(C3*exp(-Ea/(k*T1))) + (1-d)*(y*(C2*exp(-Ep/(k*T1))) + (1-y)*(C1*exp(-Ef/(k*T1))))
    N2 <- d*(C3*exp(-Ea/(k*T2))) + (1-d)*(y*(C2*exp(-Ep/(k*T2))) + (1-y)*(C1*exp(-Ef/(k*T2))))
    Q10_new <- N2/N1
    
    # percentage change from 0.65eV baseline
    perc_change <- (((Q10_new/Q10_standard)-1)*100)
    
    # work out the emergent ecosystem E
    E_eco <- d*Ea + ((1-d)*(y*Ep + (1-y)*Ef))
    
    eco_E <- c(eco_E, E_eco)
    Q10 <- c(Q10, Q10_new)
    change <- c(change, perc_change)
    
  }
  
  
}

results_df <- data.frame(auto_hetero, prok_euk, Q10, change, eco_E)

ggplot(results_df, aes(x = (1-auto_hetero)*100, y = prok_euk*100, fill = eco_E)) + 
  geom_tile() +
  scale_fill_gradient(low = "#35978f", high = "#bf812d") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  stat_contour(aes(z = Q10), bins = 4, color = "black") +
  labs(title = "Changes in Ecosystem Respiration",
       x = "Percentage Heterotrophs vs Autotrophs",
       y = "Percentage Bacteria within Heterotrophs")

short_df <- results_df[results_df$auto_hetero > 0.25 & results_df$auto_hetero < 0.75,]
short_df <- short_df[short_df$prok_euk > 0.25 & short_df$prok_euk < 0.75,]

# before plotting properly, check values for contour lines...

p <- ggplot(aes((1-auto_hetero)*100, y = prok_euk*100, z=change), data = short_df) + 
  geom_raster(data=short_df, aes(fill=change)) +
  scale_fill_gradient(limits=range(short_df$change), high = 'red', low = 'yellow')
p

# Plot 2: This plot adds the isolines but no labels and it also adds a second legend for level which I don't want
p <- p + geom_contour(aes(colour = ..level..), bins = 4, color='gray30', na.rm=T,     show.legend=T)
p

# Plot 3: This plot has the labeled isolines but it removes the z legend that I want to show
direct.label(p, list("bottom.pieces", colour='black')) 

# OK, so for all intents and purposes, contours for % change are 4, 8, 12 and 16

p1 <- ggplot(short_df, aes(x = (1-auto_hetero)*100, y = prok_euk*100)) + 
  geom_tile(aes(fill = Q10)) +
  #scale_fill_gradient(name = "Q10", low = "yellow", high = "red") +
  scale_fill_gradient(name = expression("Q"["10"]), low = "white", high = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  stat_contour(aes(z = change), bins = 4, color = "grey50", size = 2) +
  labs(title = "Short-Term Changes in Ecosystem Respiration",
       x = "Percentage Heterotrophs vs Autotrophs",
       y = "Percentage Bacteria within Heterotrophs") +
  coord_fixed() +
  main_theme

p2 <- ggplot(short_df, aes(x = (1-auto_hetero)*100, y = prok_euk*100)) + 
  geom_tile(aes(fill = change)) +
  scale_fill_gradient(name = "% Flux\nIncrease", low = "white", high = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  stat_contour(aes(z = change), bins = 4, color = "grey50", size = 2) +
  labs(title = "Short-Term Changes in Ecosystem Respiration",
       x = "Percentage Heterotrophs vs Autotrophs",
       y = "Percentage Bacteria within Heterotrophs") +
  coord_fixed() +
  main_theme

p3 <- ggplot(short_df, aes(x = (1-auto_hetero)*100, y = prok_euk*100)) + 
  geom_tile(aes(fill = eco_E)) +
  scale_fill_gradient(name = expression(paste("Emergent ", italic("E"))), low = "white", high = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  stat_contour(aes(z = change), bins = 4, color = "grey50", size = 2) +
  labs(title = "Short-Term Changes in Ecosystem Respiration",
       x = "Percentage Heterotrophs vs Autotrophs",
       y = "Percentage Bacteria within Heterotrophs") +
  coord_fixed() +
  main_theme

ggsave(file = '../../Results/figures/eco_resp_Q10.png', p1, width = 10, height = 8)
ggsave(file = '../../Results/figures/eco_resp_percent.png', p2, width = 10, height = 8)
ggsave(file = '../../Results/figures/eco_resp_E.png', p3, width = 10, height = 8)

# 3 legends Q10, E, percentage change

##########################################################################
## Set new parameter values to re-generate plot for long-term responses ##
##########################################################################

Ep <- 0.98 # average inter-specific E for mesophilic bacteria
T2 <- T1+4 # 4C change
C2 <- 2/(exp(-Ep/(k*T1))) # constant for prokaryotes

# standard Q10 with 0.65eV for comparison
Q10_standard <- (C1*exp(-Ef/(k*T2)))/(C1*exp(-Ef/(k*T1)))


Q10 <- c()
change <- c()
auto_hetero <- c()
prok_euk <- c()
eco_E <- c()


for(i in 1:100){
  
  for(j in 1:100){
    
    d <- 0.01*i # auto vs hetero
    y <- 0.01*j # prok vs euk
    
    auto_hetero <- c(auto_hetero, d)
    prok_euk <- c(prok_euk, y)
    
    # first part autotroph vs heterotroph, second part bacteria vs eukaryotes within heterotrophs
    N1 <- d*(C3*exp(-Ea/(k*T1))) + (1-d)*(y*(C2*exp(-Ep/(k*T1))) + (1-y)*(C1*exp(-Ef/(k*T1))))
    N2 <- d*(C3*exp(-Ea/(k*T2))) + (1-d)*(y*(C2*exp(-Ep/(k*T2))) + (1-y)*(C1*exp(-Ef/(k*T2))))
    Q10_new <- N2/N1
    
    # percentage change from 0.65eV baseline
    perc_change <- (((Q10_new/Q10_standard)-1)*100)
    
    # work out the emergent ecosystem E
    E_eco <- d*Ea + ((1-d)*(y*Ep + (1-y)*Ef))
    
    eco_E <- c(eco_E, E_eco)
    Q10 <- c(Q10, Q10_new)
    change <- c(change, perc_change)
    
  }
  
  
}

results_df <- data.frame(auto_hetero, prok_euk, Q10, change, eco_E)

ggplot(results_df, aes(x = (1-auto_hetero)*100, y = prok_euk*100, fill = change)) + 
  geom_tile() +
  scale_fill_gradient(low = "#35978f", high = "#bf812d") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  stat_contour(aes(z = eco_E), bins = 4, color = "black") +
  labs(title = "Changes in Ecosystem Respiration",
       x = "Percentage Heterotrophs vs Autotrophs",
       y = "Percentage Bacteria within Heterotrophs")

short_df <- results_df[results_df$auto_hetero > 0.25 & results_df$auto_hetero < 0.75,]
short_df <- short_df[short_df$prok_euk > 0.25 & short_df$prok_euk < 0.75,]

# before plotting properly, check values for contour lines...

p <- ggplot(aes((1-auto_hetero)*100, y = prok_euk*100, z=change), data = short_df) + 
  geom_raster(data=short_df, aes(fill=change)) +
  scale_fill_gradient(limits=range(short_df$change), high = 'red', low = 'yellow')
p

# Plot 2: This plot adds the isolines but no labels and it also adds a second legend for level which I don't want
p <- p + geom_contour(aes(colour = ..level..), bins = 3, color='gray30', na.rm=T,     show.legend=T)
p

# Plot 3: This plot has the labeled isolines but it removes the z legend that I want to show
direct.label(p, list("bottom.pieces", colour='black')) 


p4 <- ggplot(short_df, aes(x = (1-auto_hetero)*100, y = prok_euk*100)) + 
  geom_tile(aes(fill = change)) +
  scale_fill_gradient(name = "% Flux\nIncrease", low = "white", high = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  stat_contour(aes(z = change), bins = 3, color = "grey50", size = 2) +
  labs(title = "Long-Term Changes in Ecosystem Respiration",
       x = "Percentage Heterotrophs vs Autotrophs",
       y = "Percentage Bacteria within Heterotrophs") +
  coord_fixed() +
  main_theme

ggsave(file = '../../Results/figures/long_term_change.png', p4, width = 10, height = 8)

###################################################################################


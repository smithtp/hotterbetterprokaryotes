# Want to look at whether E from growth rates is similar to E from metabolic fluxes
# are there any pairwise comparisons in the data?

rm(list = ls())

require(ggplot2)
require(reshape)
require(grid)
require(gridExtra)

main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        plot.title = element_text(size=18, vjust=1),
        legend.text=element_text(size=18),
        legend.title = element_text(size = 22),
        strip.text.x = element_text(size = 16))


setwd("~/Documents/hotterbetterprokaryotes/Code/R_Scripts/")

##################################################################################
#### first we'll just look at the E from flux data against E from growth data ####
##################################################################################

# this is the flux data only
group_means <- read.csv("../../Data/summaries/group_means_non_aggregated_fluxes.csv")
group_means <- group_means[group_means$ConKingdom %in% c('Archaea', 'Bacteria'),]

# this is the growth data
data_agg <- read.csv('../../Data/summaries/all_means_non_aggregated.csv')
data_agg <- data_agg[data_agg$ConKingdom %in% c('Archaea', 'Bacteria'),]


dodge = position_dodge(width=0.3)

# we'll put the Es from the flux data in the background for this plot
Flux_back <- group_means[, c('Species', 'ConKingdom', 'ConKingdom' , 'E')]
Growth_back <- data_agg[, c('Species', 'ConKingdom', 'ConKingdom' , 'E')]
Kingdom_fluxES <- unique(group_means[, c('ConKingdom', 'ConKingdom', 'C_ES', 'C_ES_min', 'C_ES_max')])
Kingdom_growthES <- unique(data_agg[, c('ConKingdom', 'ConKingdom', 'ConKingdom_ES', 'ConKingdom_ES_min', 'ConKingdom_ES_max')])
Kingdom_fluxES$var <- 'Flux ES'
Kingdom_growthES$var <- 'Growth ES'

colnames(Kingdom_fluxES) <- c('ConKingdom', 'ConKingdom', 'E', 'E_min', 'E_max', 'var')
colnames(Kingdom_growthES) <- c('ConKingdom', 'ConKingdom', 'E', 'E_min', 'E_max', 'var')

Kingdom_front <- rbind(Kingdom_fluxES, Kingdom_growthES)

Kingdom_Overlap_Plot <- ggplot(Kingdom_front, aes(x=ConKingdom)) +
  geom_errorbar(aes(ymax=E_max, ymin=E_min, group=var), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E, colour=var, group=var, shape = var), size = 4, position=dodge) +  
  ylab('Thermal Sensitivity (E, eV)') +
  xlab('Kingdom') +
  geom_hline(aes(yintercept=0.65), linetype='dotted') +
  ylim(0, 3) +
  main_theme +
  scale_colour_manual(guide='legend', values=c("#CC0066", "#339900"),
                      labels = c(expression(italic('Flux E'['S'])),
                                 expression(italic('Growth E'['S'])))) +
  scale_shape_discrete(guide='legend',
                        labels = c(expression(italic('Flux E'['S'])),
                                   expression(italic('Growth E'['S'])))) +
  theme(legend.title=element_blank(),
        legend.justification=c(0,2), legend.position=c(0.1,1),
        plot.margin = unit(c(1,0,1,1), "cm")) +
  annotate("text", x = 0.5, y = 3, label = "A", size = 10)
Kingdom_Overlap_Plot


## lets add a density plot alongside

Flux_Density_Plot <- ggplot(Flux_back, aes(x=E)) +
  geom_density(aes(group=ConKingdom, fill=ConKingdom), alpha=0.5) +
  xlim(0, 3) +
  ylab("Density") +
  geom_vline(aes(xintercept=0.65), linetype='dotted') +
  coord_flip() + scale_y_reverse() +
  main_theme +
  scale_fill_manual(guide='legend',
                    name = "Flux",
                      values =c('Archaea'='#FF33FF',
                                'Bacteria'='#CC0066'),
                      labels = c("Archaea",
                                 "Bacteria")) +
  theme(legend.justification=c(0,2), legend.position=c(0.1,1.1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin = unit(c(1,0,1,0), "cm")) +
  annotate("text", x = 3, y = 1.9, label = "B", size = 10)
Flux_Density_Plot

Growth_Density_Plot <- ggplot(Growth_back, aes(x=E)) +
  geom_density(aes(group=ConKingdom, fill=ConKingdom), alpha=0.5) +
  ylab("Density") +
  xlim(0, 3) +
  geom_vline(aes(xintercept=0.65), linetype='dotted') +
  coord_flip() + scale_y_reverse() +
  main_theme +
  scale_fill_manual(guide='legend',
                    name = "Growth",
                    values =c('Archaea'='#66CC00',
                              'Bacteria'='#009933'),
                    labels = c("Archaea",
                               "Bacteria")) +
  theme(legend.justification=c(0,2), legend.position=c(0.1,1.1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin = unit(c(1,0,1,0), "cm")) +
  annotate("text", x = 3, y = 1, label = "C", size = 10)
Growth_Density_Plot

# now output the plot
png('../../Results/figures/Fluxes_versus_Growth_E.png', width = 960, height = 480)
grid.arrange(Kingdom_Overlap_Plot, Flux_Density_Plot, Growth_Density_Plot, ncol=3, nrow=1, widths=c(2, 2, 2), heights=4, clip = TRUE)
dev.off()



### Add in data from Sofia on plant resp fluxes

e_data <- read.csv("../../Data/sofia_data/DataAutotrophs_nopseudo.csv")

# density plot

sofia_density_plot <- ggplot(e_data, aes(x = E)) +
  geom_density(aes(group = id_process, fill = id_process), alpha = 0.5) +
  xlim(0, 3) +
  ylab("Density") +
  geom_vline(aes(xintercept=0.65), linetype='dotted') +
  coord_flip() + scale_y_reverse() +
  main_theme +
  scale_fill_manual(guide='legend',
                    name = "Autotroph Flux",
                    values =c('respiration rate'='#CC3300'),
                    labels = c("Respiration")) +
  theme(legend.justification=c(0,2), legend.position=c(0.1,1.05),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin = unit(c(1,1,1,0), "cm")) +
  annotate("text", x = 3, y = 1.5, label = "D", size = 10)
sofia_density_plot


## now add this to our previous plot panel as a nice comparison

png('../../Results/figures/Fluxes_versus_Growth_E_extra.png', width = 1280, height = 480)
grid.arrange(Kingdom_Overlap_Plot, Flux_Density_Plot, Growth_Density_Plot, sofia_density_plot, ncol=4, nrow=1, widths=c(2, 2, 2, 2), heights=4, clip = TRUE)
dev.off()

# and just check what the mean/median of this data is:

mean(e_data$E_sch)

median(e_data$E_sch)

# implement bootstrapping code for CIs, like the python code

resample_all <- function(x){
  boot_nums = c()
  len_x = length(x)
  for (i in 1:5000){
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

bootstrap_upper(e_data$E_sch)
bootstrap_mean(e_data$E_sch)
bootstrap_lower(e_data$E_sch)


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
data_agg <- read.csv('../../Data/summaries/all_means_non_aggregated.csv')
data_agg <- data_agg[data_agg$ConKingdom %in% c('Archaea', 'Bacteria'),]

#Sys.setlocale("LC_CTYPE", "Latvian") #Needed to make E appear (in windows?)
dodge = position_dodge(width=0.18)

Kingdom_data <- data_agg[data_agg$ConKingdom_len > 4, ]

Kingdom_back <- Kingdom_data[, c('Species', 'ConKingdom', 'ConKingdom' , 'E')]
Kingdom_ES <- unique(Kingdom_data[, c('ConKingdom', 'ConKingdom', 'ConKingdom_ES', 'ConKingdom_ES_min', 'ConKingdom_ES_max')])
Kingdom_EG <- unique(Kingdom_data[, c('ConKingdom', 'ConKingdom', 'ConKingdom_EG', 'ConKingdom_EG_min', 'ConKingdom_EG_max')])
Kingdom_ES$var <- 'ES'
Kingdom_EG$var <- 'EG'

colnames(Kingdom_ES) <- c('ConKingdom', 'ConKingdom', 'E', 'E_min', 'E_max', 'var')
colnames(Kingdom_EG) <- c('ConKingdom', 'ConKingdom', 'E', 'E_min', 'E_max', 'var')

Kingdom_front <- rbind(Kingdom_ES, Kingdom_EG)

Kingdom_Overlap_Plot <- ggplot(Kingdom_front, aes(x=ConKingdom)) +
  geom_point(data=Kingdom_back, aes(y=E, colour='3', shape = '3'), alpha = 0.7, position = position_jitter(width = 0.4, height = 0.0)) +
  geom_errorbar(aes(ymax=E_max, ymin=E_min, group=var), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E, colour=var, shape = var, group=var), size = 4, position=dodge) +  
  xlab('Kingdom') +
  ylab('Activation Energy (E)') +
  geom_hline(aes(yintercept=0.65), linetype='dotted') +
  scale_colour_manual(guide='legend',
                      values =c('3'='#E69F00',
                                'ES'='#D55E00',
                                'EG'='royalblue3'),
                      labels = c(expression(italic('E'['S'])),
                                 expression(italic('E'['G'])),
                                 quote(italic('\U0112'[S])))
  ) +
  scale_shape_manual(guide = 'legend',
                     values =c('3'=19,
                               'ES'=15,
                               'EG'=17),
                     labels = c(expression(italic('E'['S'])),
                                expression(italic('E'['G'])),
                                quote(italic('\U0112'[S])))
  ) +
  coord_cartesian(ylim = c(0,5)) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.justification=c(0,1), legend.position=c(0,1),
        plot.margin = unit(c(1,0,1,1), "cm"))
Kingdom_Overlap_Plot

## lets try to add a density plot alongside

Kingdom_Density_Plot <- ggplot(Kingdom_back, aes(x=E)) +
  geom_density(aes(group=ConKingdom, colour=ConKingdom, fill=ConKingdom), alpha=0.3) +
  ylab("Density") +
  geom_vline(aes(xintercept=0.65), linetype='dotted') +
  coord_flip(xlim = c(0, 5)) + scale_y_reverse() +
  main_theme +
  theme(legend.title=element_blank(),
        legend.justification=c(0,1), legend.position=c(0,1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin = unit(c(1,1,1,0), "cm"))
Kingdom_Density_Plot

## TempPref stuff here

data_agg <- data_agg[!is.na(data_agg$TempPref),]
TempPref_back <- data_agg[, c('Species', 'ConKingdom', 'TempPref' , 'E')]

TempPref_back$TempPref <- as.character(TempPref_back$TempPref)
#TempPref_back$TempPref <- factor(TempPref_back$TempPref, levels=c("Thermophile", "Mesophile"))

TempPref_ES <- unique(data_agg[, c('ConKingdom', 'TempPref', 'TempPref_ES', 'TempPref_ES_min', 'TempPref_ES_max')])
TempPref_EG <- unique(data_agg[, c('ConKingdom', 'TempPref', 'TempPref_EG', 'TempPref_EG_min', 'TempPref_EG_max')])
TempPref_ES$var <- 'ES'
TempPref_EG$var <- 'EG'

colnames(TempPref_ES) <- c('ConKingdom', 'TempPref', 'E', 'E_min', 'E_max', 'var')
colnames(TempPref_EG) <- c('ConKingdom', 'TempPref', 'E', 'E_min', 'E_max', 'var')

TempPref_front <- rbind(TempPref_ES, TempPref_EG)

TempPref_front$TempPref <- as.character(TempPref_front$TempPref)
#TempPref_front$TempPref <- factor(TempPref_front$TempPref, levels=c("Thermophile", "Mesophile"))

TempPref_Overlap_Plot <- ggplot(TempPref_front, aes(x=TempPref)) +
  geom_point(data=TempPref_back, aes(y=E, colour='3', shape ='3'), alpha = 0.7, position = position_jitter(width = 0.3, height = 0.0)) +
  geom_errorbar(aes(ymax=E_max, ymin=E_min, group=var), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E, colour=var, group=var, shape=var), size=4, position=dodge) +  
  ylab('Thermal Sensitivity (E, eV)') +
  #xlab('Temperature Preference') +
  geom_hline(aes(yintercept=0)) +
  #geom_hline(aes(yintercept=0.65), linetype='dotted') +
  scale_colour_manual(guide=FALSE,
                      values =c('3'='#E69F00',
                                'ES'='#D55E00',
                                'EG'='royalblue3')
  ) +
  scale_shape_manual(guide = FALSE,
                     values =c('3'=19,
                               'ES'=15,
                               'EG'=17)
  ) +
  facet_wrap(~ConKingdom, ncol=2, scales='free') +
  coord_cartesian(ylim=c(-1, 5)) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.position='bottom',
        axis.title.x = element_blank()) 
TempPref_Overlap_Plot


ggsave(file = '../../Results/without_aggregation/BA_Fits/1_Kingdom_Overlap_Plot.png', Kingdom_Overlap_Plot, width = 12, height = 6)

png('../../Results/figures/Kingdom_Overlap_Plot_extra.png', width = 960, height = 480)
grid.arrange(Kingdom_Overlap_Plot, Kingdom_Density_Plot, ncol=2, nrow=1, widths=c(4, 1), heights=4, clip = TRUE)
dev.off()

ggsave(file = '../../Results/figures/1_Temp_Overlap_Plot_2.png', TempPref_Overlap_Plot, width = 16, height = 4)

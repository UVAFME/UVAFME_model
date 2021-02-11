#
# #########################
# Purpose: Test UVAFME output 
# Author: Adrianna C. Foster
# Date: October, 2018
# R version 3.5.1 (2018-07-02) 'Feather Spray'
# #########################
# #########################
# Import format: .csv
# #########################

library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)

## Output directory
output_dir <- 'output_data_sample'

## Live species-level output 
spec_dat <- read.csv(paste0(output_dir, '/Species_Data.csv'), 
                     stringsAsFactors = FALSE) 
# Species not present at site coded as -999
spec_dat[spec_dat==-999] <- NA

## Live genus-level output 
gen_dat <- read.csv(paste0(output_dir, '/Genus_Data.csv'), 
                     stringsAsFactors = FALSE) 
# Genera not present at site coded as -999
gen_dat[gen_dat==-999] <- NA

## Dead species-level output 
dead_spec_dat <- read.csv(paste0(output_dir, '/Dead_Species_Data.csv'), 
                     stringsAsFactors = FALSE) 
# Species not present at site coded as -999
dead_spec_dat[dead_spec_dat==-999] <- NA

## Dead genus-level output 
dead_gen_dat <- read.csv(paste0(output_dir, '/Dead_Genus_Data.csv'), 
                    stringsAsFactors = FALSE) 
# Genera not present at site coded as -999
dead_gen_dat[dead_gen_dat==-999] <- NA

## Across-species output 
tot_dat <- read.csv(paste0(output_dir, '/Total_Plot_Values.csv'), 
                     stringsAsFactors = FALSE) 

## Climate/site output 
clim_dat <- read.csv(paste0(output_dir, '/Climate.csv'), 
                    stringsAsFactors = FALSE) 

## Soil output 
soil_dat <- read.csv(paste0(output_dir, '/SoilDecomp.csv'), 
                     stringsAsFactors = FALSE) 

## Using ggplot and patchwork here 
ggsp <- ggplot(spec_dat, aes(year, total_biomC)) + 
  geom_area(aes(fill = species), colour = 'white') + 
  facet_wrap(.~siteID) + 
  scale_fill_manual(values = c('lightsalmon', 'lightseagreen',
                               'navyblue', 'gold'),
                    labels = c('Alaska birch', 'White spruce',
                               'Black spruce', 'Quaking aspen'),
                    name = 'Species') +
  xlab('Stand Age (years)') + 
  ylab(expression(Aboveground~Biomass~(tC~ha^-1))) +
  theme_bw()

gggen <- ggplot(gen_dat, aes(year, total_biomC)) + 
  geom_area(aes(fill = genus), colour = 'white') + 
  facet_wrap(.~siteID) + 
  scale_fill_manual(values = c('lightsalmon', 'seagreen4', 
                               'gold'),
                    labels = c('Birch', 'Spruce', 'Aspen'),
                    name = 'Genus') +
  xlab('Stand Age (years)') + 
  ylab(expression(Aboveground~Biomass~(tC~ha^-1))) +
  theme_bw()

gdsp <- ggplot(dead_spec_dat, aes(year, total_biomC)) + 
  geom_area(aes(fill = species), colour = 'white') + 
  facet_wrap(.~siteID) + 
  scale_fill_manual(values = c('lightsalmon', 'lightseagreen',
                               'navyblue', 'gold'),
                    labels = c('Alaska birch', 'White spruce',
                               'Black spruce', 'Quaking aspen'),
                    name = 'Species') +
  xlab('Stand Age (years)') + 
  ylab(expression(Dead~Aboveground~Biomass~(tC~ha^-1))) +
  theme_bw()

ggdgen <- ggplot(dead_gen_dat, aes(year, total_biomC)) + 
  geom_area(aes(fill = genus), colour = 'white') + 
  facet_wrap(.~siteID) + 
  scale_fill_manual(values = c('lightsalmon', 'seagreen4', 
                               'gold'),
                    labels = c('Birch', 'Spruce', 'Aspen'),
                    name = 'Genus') +
  xlab('Stand Age (years)') + 
  ylab(expression(Dead~Aboveground~Biomass~(tC~ha^-1))) +
  theme_bw()

ggtot <- ggplot(tot_dat, aes(year, total_biomC)) +
  geom_ribbon(aes(x = year, ymin = total_biomC - total_biomC_sd,
                  ymax = total_biomC + total_biomC_sd,
                  fill = as.factor(siteID)),
              alpha = 0.5) + 
  geom_line(aes(colour = as.factor(siteID)), size = 1) +
  scale_colour_discrete(name = 'siteID') +
  scale_fill_discrete(name = 'siteID') +
  xlab('Stand Age (years)') + 
  ylab(expression(Total~Biomass~"±"~SD~(tC~ha^-1))) + 
  theme_bw()

ggclim <- ggplot(clim_dat, aes(year, pet)) +
  geom_line(aes(colour = as.factor(siteID)), size = 1, show.legend = FALSE) +
  scale_colour_discrete(name = 'siteID') +
  xlab('Stand Age (years)') + 
  ylab('Annual Potential Evapotranspiration (cm)') + 
  theme_bw()


ggorg <- ggplot(soil_dat, aes(year, odepth)) +
  geom_ribbon(aes(x = year, ymin = odepth - odepth_sd,
                  ymax = odepth + odepth_sd,
                  fill = as.factor(siteID)),
              alpha = 0.5) + 
  geom_line(aes(colour = as.factor(siteID)), size = 1) +
  scale_colour_discrete(name = 'siteID') +
  scale_fill_discrete(name = 'siteID') +
  xlab('Stand Age (years)') + 
  ylab(expression(Organic~Layer~Depth~"±"~SD~(cm))) + 
  theme_bw()

p1 <- ggsp + gggen + gdsp + ggtot +
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(title = 'Vegetation Output from UVAFME')

p2 <- ggclim + ggorg + 
  plot_annotation(title = 'Site Condition Output from UVAFME')


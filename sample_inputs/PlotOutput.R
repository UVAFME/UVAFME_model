#
# #########################
# Purpose: Plot some output
# Author: Adrianna C. Foster
# Date: August, 2022
# R version 4.1.1 (2022-04-22) 'Vigorous Calisthenics'
# #########################
# #########################
# Input format: csv
# #########################

library(dplyr)
library(ggplot2)

setwd('/Users/afoster/Documents/ABoVE/UVAFME_model/')

specDat <- read.csv('sample_inputs/output_data_sample/Species_Data.csv')
specDat[specDat==-999] <- NA

clim <- read.csv('sample_inputs/output_data_sample/Climate.csv')

ggplot(specDat, aes(year, total_biomC)) +
  geom_area(aes(fill = species)) +
  facet_wrap(.~siteID)

ggplot(clim, aes(year, organic_depth)) +
  geom_line(aes(colour = as.factor(siteID)))

#
# #########################
# Purpose: Create climate input data for UVAFME
# Author: Adrianna C. Foster
# Date: August, 2022
# R version 4.1.1 (2022-04-22) 'Vigorous Calisthenics'
# #########################
# #########################
# Input format: csv
# Output format: csv
# #########################

library(dplyr)
library(reshape2)
library(rgdal)
library(raster)
library(data.table)

setwd('/Users/afoster/Documents/ABoVE/UVAFME_model/')


### Functions ------------------------------------------------------------------

parameterize_climate <- function(climNA_fname, cldmn_dir, cldsd_dir,
                                 wind_dir){
  ## function to parameterize climate files from ClimateNA output and wind/cld
  ## maps
  ## Inputs:
  ##  climNA_fname - monthly historical time series output from ClimateNA
  ##  cldmn_dir - directory with mean monthly cloudiness values
  ##  cldsd_dir - directory with sd of monthly cloudiness values
  ##  wind_dir - directory with mean monthly wind speed
  ## Outputs:
  ##  climMn  - data frame for _climate.csv file
  ##  climSd  - data frame for _climate_stddev.csv file 
  ##  extraMN - data frame for _climate_ex.csv file
  ##  extraSd - data frame for _climate_ex_stddev.csv file
  
  climNA <- fread(climNA_fname)
  
  #get rid of columns we don't want
  climNA <- dplyr::select(climNA, Year:Longitude, Tmax01:Tmin12, PPT01:PPT12, 
                          RH01:RH12)
  
  
  #summarize data - mean
  climNA_means <- dplyr::group_by(climNA, ID2) %>%
    dplyr::summarize(tmin_jan = mean(Tmin01, na.rm = TRUE), 
                     tmin_feb = mean(Tmin02, na.rm = TRUE),
                     tmin_mar = mean(Tmin03, na.rm = TRUE), 
                     tmin_apr = mean(Tmin04, na.rm = TRUE),
                     tmin_may = mean(Tmin05, na.rm = TRUE), 
                     tmin_jun = mean(Tmin06, na.rm = TRUE), 
                     tmin_jul = mean(Tmin07, na.rm = TRUE), 
                     tmin_aug = mean(Tmin08, na.rm = TRUE),
                     tmin_sep = mean(Tmin09, na.rm = TRUE),
                     tmin_oct = mean(Tmin10, na.rm = TRUE),
                     tmin_nov = mean(Tmin11, na.rm = TRUE),
                     tmin_dec = mean(Tmin12, na.rm = TRUE),
                     
                     tmax_jan = mean(Tmax01, na.rm = TRUE), 
                     tmax_feb = mean(Tmax02, na.rm = TRUE),
                     tmax_mar = mean(Tmax03, na.rm = TRUE), 
                     tmax_apr = mean(Tmax04, na.rm = TRUE),
                     tmax_may = mean(Tmax05, na.rm = TRUE), 
                     tmax_jun = mean(Tmax06, na.rm = TRUE), 
                     tmax_jul = mean(Tmax07, na.rm = TRUE), 
                     tmax_aug = mean(Tmax08, na.rm = TRUE),
                     tmax_sep = mean(Tmax09, na.rm = TRUE),
                     tmax_oct = mean(Tmax10, na.rm = TRUE),
                     tmax_nov = mean(Tmax11, na.rm = TRUE),
                     tmax_dec = mean(Tmax12, na.rm = TRUE),
                     
                     prcp_jan = mean(PPT01, na.rm = TRUE), 
                     prcp_feb = mean(PPT02, na.rm = TRUE),
                     prcp_mar = mean(PPT03, na.rm = TRUE), 
                     prcp_apr = mean(PPT04, na.rm = TRUE),
                     prcp_may = mean(PPT05, na.rm = TRUE), 
                     prcp_jun = mean(PPT06, na.rm = TRUE), 
                     prcp_jul = mean(PPT07, na.rm = TRUE), 
                     prcp_aug = mean(PPT08, na.rm = TRUE),
                     prcp_sep = mean(PPT09, na.rm = TRUE),
                     prcp_oct = mean(PPT10, na.rm = TRUE),
                     prcp_nov = mean(PPT11, na.rm = TRUE),
                     prcp_dec = mean(PPT12, na.rm = TRUE),
                     
                     rh_jan = mean(RH01, na.rm = TRUE), 
                     rh_feb = mean(RH02, na.rm = TRUE),
                     rh_mar = mean(RH03, na.rm = TRUE), 
                     rh_apr = mean(RH04, na.rm = TRUE),
                     rh_may = mean(RH05, na.rm = TRUE), 
                     rh_jun = mean(RH06, na.rm = TRUE), 
                     rh_jul = mean(RH07, na.rm = TRUE), 
                     rh_aug = mean(RH08, na.rm = TRUE),
                     rh_sep = mean(RH09, na.rm = TRUE),
                     rh_oct = mean(RH10, na.rm = TRUE),
                     rh_nov = mean(RH11, na.rm = TRUE),
                     rh_dec = mean(RH12, na.rm = TRUE))
  
  
  #summarize data - sd
  climNA_sds <- group_by(climNA, ID2) %>%
    dplyr::summarize(tmin_std_jan = sd(Tmin01, na.rm = TRUE), 
                     tmin_std_feb = sd(Tmin02, na.rm = TRUE),
                     tmin_std_mar = sd(Tmin03, na.rm = TRUE), 
                     tmin_std_apr = sd(Tmin04, na.rm = TRUE),
                     tmin_std_may = sd(Tmin05, na.rm = TRUE), 
                     tmin_std_jun = sd(Tmin06, na.rm = TRUE), 
                     tmin_std_jul = sd(Tmin07, na.rm = TRUE), 
                     tmin_std_aug = sd(Tmin08, na.rm = TRUE),
                     tmin_std_sep = sd(Tmin09, na.rm = TRUE),
                     tmin_std_oct = sd(Tmin10, na.rm = TRUE),
                     tmin_std_nov = sd(Tmin11, na.rm = TRUE),
                     tmin_std_dec = sd(Tmin12, na.rm = TRUE),
                     
                     tmax_std_jan = sd(Tmax01, na.rm = TRUE), 
                     tmax_std_feb = sd(Tmax02, na.rm = TRUE),
                     tmax_std_mar = sd(Tmax03, na.rm = TRUE), 
                     tmax_std_apr = sd(Tmax04, na.rm = TRUE),
                     tmax_std_may = sd(Tmax05, na.rm = TRUE), 
                     tmax_std_jun = sd(Tmax06, na.rm = TRUE), 
                     tmax_std_jul = sd(Tmax07, na.rm = TRUE), 
                     tmax_std_aug = sd(Tmax08, na.rm = TRUE),
                     tmax_std_sep = sd(Tmax09, na.rm = TRUE),
                     tmax_std_oct = sd(Tmax10, na.rm = TRUE),
                     tmax_std_nov = sd(Tmax11, na.rm = TRUE),
                     tmax_std_dec = sd(Tmax12, na.rm = TRUE),
                     
                     prcp_std_jan = sd(PPT01, na.rm = TRUE), 
                     prcp_std_feb = sd(PPT02, na.rm = TRUE),
                     prcp_std_mar = sd(PPT03, na.rm = TRUE), 
                     prcp_std_apr = sd(PPT04, na.rm = TRUE),
                     prcp_std_may = sd(PPT05, na.rm = TRUE), 
                     prcp_std_jun = sd(PPT06, na.rm = TRUE), 
                     prcp_std_jul = sd(PPT07, na.rm = TRUE), 
                     prcp_std_aug = sd(PPT08, na.rm = TRUE),
                     prcp_std_sep = sd(PPT09, na.rm = TRUE),
                     prcp_std_oct = sd(PPT10, na.rm = TRUE),
                     prcp_std_nov = sd(PPT11, na.rm = TRUE),
                     prcp_std_dec = sd(PPT12, na.rm = TRUE),
                     
                     rh_std_jan = sd(RH01, na.rm = TRUE), 
                     rh_std_feb = sd(RH02, na.rm = TRUE),
                     rh_std_mar = sd(RH03, na.rm = TRUE), 
                     rh_std_apr = sd(RH04, na.rm = TRUE),
                     rh_std_may = sd(RH05, na.rm = TRUE), 
                     rh_std_jun = sd(RH06, na.rm = TRUE), 
                     rh_std_jul = sd(RH07, na.rm = TRUE), 
                     rh_std_aug = sd(RH08, na.rm = TRUE),
                     rh_std_sep = sd(RH09, na.rm = TRUE),
                     rh_std_oct = sd(RH10, na.rm = TRUE),
                     rh_std_nov = sd(RH11, na.rm = TRUE),
                     rh_std_dec = sd(RH12, na.rm = TRUE))  
  
  #get site lat and long
  site_locs <- dplyr::select(climNA, ID2, Latitude, Longitude) %>% 
    distinct_all()
  
  #merge climate data with site locations to get latitude/longitude
  climNA_means <- merge(climNA_means, site_locs, by = 'ID2')
  climNA_sds <- merge(climNA_sds, site_locs, by = 'ID2')
  
  #rename columns and reorder
  climNA_means <- dplyr::rename(climNA_means, 'site' = 'ID2', 'latitude' = 'Latitude',
                                'longitude' = 'Longitude') %>% 
    dplyr::select(site, latitude, longitude, tmin_jan:rh_dec)
  
  climNA_sds <- dplyr::rename(climNA_sds, 'site' = 'ID2', 'latitude' = 'Latitude',
                              'longitude' = 'Longitude') %>% 
    dplyr::select(site, latitude, longitude, tmin_std_jan:rh_std_dec)
  
  #Monthly cloudiness input---
  
  #get input cloud tif files
  cldmns <- list.files(cldmn_dir, pattern = '.tif$')
  cldsds <- list.files(cldsd_dir, pattern = '.tif$')
  
  #read in stack of cloud means and sds
  cldmn <- stack()
  for (r in 1:length(cldmns)){
    cld_tmp <- raster(paste(cldmn_dir, cldmns[r], sep = '/'))
    cldmn <- stack(cldmn, cld_tmp)
  }
  
  cldsd <- stack()
  for (r in 1:length(cldsds)){
    cld_tmp <- raster(paste(cldsd_dir, cldsds[r], sep = '/'))
    cldsd <- stack(cldsd, cld_tmp)
  }
  
  
  #convert site info to spatial points dataframe
  siteLocs_sp <- SpatialPointsDataFrame(coords = site_locs[,c(3,2)], 
                                        data = site_locs, 
                                        proj4string = crs(cldmn))
  
  #extract cloud information at site locations
  ex_mn <- extract(cldmn, siteLocs_sp)
  ex_sd <- extract(cldsd, siteLocs_sp)
  siteLocs_sp@data <- cbind(siteLocs_sp@data, data.frame(ex_mn), 
                            data.frame(ex_sd))
  
  #convert back to dataframe and split into mean and std columns
  cld_mn <- data.frame(siteLocs_sp) %>% dplyr::select(ID2:cld_12mean)
  cld_sd <- data.frame(siteLocs_sp) %>% dplyr::select(ID2:Longitude, 
                                                      cld_01std:cld_12std)
  
  #rename columns
  cld_mn <- dplyr::rename(cld_mn, 'site' = 'ID2', 'latitude' = 'Latitude', 
                          'longitude' = 'Longitude', 'cld_jan' = 'cld_01mean',
                          'cld_feb' = 'cld_02mean', 'cld_mar' = 'cld_03mean',
                          'cld_apr' = 'cld_04mean', 'cld_may' = 'cld_05mean',
                          'cld_jun' = 'cld_06mean', 'cld_jul' = 'cld_07mean',
                          'cld_aug' = 'cld_08mean', 'cld_sep' = 'cld_09mean',
                          'cld_oct' = 'cld_10mean', 'cld_nov' = 'cld_11mean',
                          'cld_dec' = 'cld_12mean')
  
  cld_sd <- dplyr::rename(cld_sd, 'site' = 'ID2', 'latitude' = 'Latitude', 
                          'longitude' = 'Longitude', 'cld_jan_std' = 'cld_01std',
                          'cld_feb_std' = 'cld_02std', 'cld_mar_std' = 'cld_03std',
                          'cld_apr_std' = 'cld_04std', 'cld_may_std' = 'cld_05std',
                          'cld_jun_std' = 'cld_06std', 'cld_jul_std' = 'cld_07std',
                          'cld_aug_std' = 'cld_08std', 'cld_sep_std' = 'cld_09std',
                          'cld_oct_std' = 'cld_10std', 'cld_nov_std' = 'cld_11std',
                          'cld_dec_std' = 'cld_12std')
  
  
  #Monthly wind speed input---
  
  #get input wind tif files
  windmns <- list.files(wind_dir, pattern = '.tif$')
  
  #read in stack of wind speed means
  windmn <- stack()
  for (r in 1:length(windmns)){
    wind_tmp <- raster(paste(wind_dir, windmns[r], sep = '/'))
    windmn <- stack(windmn, wind_tmp)
  }
  
  #convert site info to spatial points dataframe
  siteLocs_sp <- SpatialPointsDataFrame(coords = site_locs[,c(3,2)], 
                                        data = site_locs, 
                                        proj4string = crs(windmn))
  
  #extract cloud information at site locations
  ex_mn <- extract(windmn, siteLocs_sp)
  siteLocs_sp@data <- cbind(siteLocs_sp@data, data.frame(ex_mn))
  
  #convert back to dataframe and split into mean and std columns
  wind_mn <- data.frame(siteLocs_sp) %>% dplyr::select(ID2:wc2.1_30s_wind_12)
  
  #rename columns
  wind_mn <- dplyr::rename(wind_mn, 'site' = 'ID2', 'latitude' = 'Latitude', 
                           'longitude' = 'Longitude', 'wind_jan' = 'wc2.1_30s_wind_01',
                           'wind_feb' = 'wc2.1_30s_wind_02', 'wind_mar' = 'wc2.1_30s_wind_03',
                           'wind_apr' = 'wc2.1_30s_wind_04', 'wind_may' = 'wc2.1_30s_wind_05',
                           'wind_jun' = 'wc2.1_30s_wind_06', 'wind_jul' = 'wc2.1_30s_wind_07',
                           'wind_aug' = 'wc2.1_30s_wind_08', 'wind_sep' = 'wc2.1_30s_wind_09',
                           'wind_oct' = 'wc2.1_30s_wind_10', 'wind_nov' = 'wc2.1_30s_wind_11',
                           'wind_dec' = 'wc2.1_30s_wind_12')
  
  extraMN <- merge(cld_mn, wind_mn, by = c('site', 'latitude', 'longitude'))
  extraMN <- merge(extraMN, dplyr::select(climNA_means, site:longitude, 
                                          rh_jan:rh_dec), 
                   by = c('site', 'latitude', 'longitude'))
  extraMN <- dplyr::select(extraMN, site:longitude, cld_jan:cld_dec, 
                           rh_jan:rh_dec, wind_jan:wind_dec)
  
  extraSD <- merge(cld_sd, dplyr::select(climNA_sds, site:longitude, 
                                         rh_std_jan:rh_std_dec), 
                   by = c('site', 'latitude', 'longitude'))
  
  climNA_means <- dplyr::select(climNA_means, site:prcp_dec)
  climNA_sds <- dplyr::select(climNA_sds, site:prcp_std_dec)
  
  
  return(list(climMn = climNA_means, climSd = climNA_sds, extraMN = extraMN,
              extraSd = extraSD))
}
parameterize_lightning <- function(lightning_dir, lightning_file, locs){
  ## function to parameterize lightning files from map
  ## Inputs:
  ##  lightning_dir  - directory with mean monthly lightning strike density
  ##  lightning_file - file with mean monthly lightning strike density
  ##  locs           - site locations
  ## Outputs:
  ##  lightning_dat  - data frame for _lightning.csv file
  
  shp <- readOGR(dsn = lightning_dir, layer = lightning_file)
  
  #convert site locations to a spatial points dataframe
  siteLocs_sp <- SpatialPointsDataFrame(coords = 
                                          locs[,c('longitude','latitude')], 
                                        data = locs,
                                        proj4string = CRS("+init=epsg:4326"))
  
  #extract site location data and convert NAs (i.e. data not present) to 0
  ex <- over(spTransform(siteLocs_sp, crs(shp)), shp)
  
  #cbind with original data frame
  siteLocs_sp@data <- cbind(siteLocs_sp@data, ex)
  
  lightning_dat <- as.data.frame(siteLocs_sp)
  
  ## get rid of unwanted columns, set NA values to 0
  lightning_dat <- dplyr::select(lightning_dat, site:longitude, 
                                 strmn_jan:strsd_dec)
  
  lightning_dat[is.na(lightning_dat)] <- 0.0
  
  
  return(lightning_dat)
  
}

## output from Climate NA 1960 - 1990 Historical Series - Monthly Variables
climNA_fname <- '../climate_dat/sample_outputs_CLIMNA_1960-1990M.csv'

## cloud cover
cldmn_dir <- paste0('/Users/afoster/Documents/ABoVE/climate_dat/cld')
cldsd_dir <- paste0('/Users/afoster/Documents/ABoVE/climate_dat/cld_std')

## wind speed
wind_dir <- paste0('/Users/afoster/Documents/ABoVE/climate_dat/wind/',
                   'wc2.1_30s_wind/')

## lightning
lightning_dir <- paste0('/Users/afoster/Documents/ABoVE/climate_dat/lightning/')
lightning_file <- 'AKLDN_10km'

## create UVAFME input data frames
clim_files <- parameterize_climate(climNA_fname = climNA_fname,
                                   cldmn_dir = cldmn_dir,
                                   cldsd_dir = cldsd_dir,
                                   wind_dir = wind_dir)

locs <- dplyr::select(clim_files$climMn, site, latitude, longitude)

lightning_dat <- parameterize_lightning(lightning_dir, lightning_file, locs)

## write to file
out.dir <- 'sample_inputs/input_data_sample'

write.csv(clim_files$climMn, paste0(out.dir, '/UVAFME2018_climate.csv'),
          row.names = FALSE, quote = FALSE)
write.csv(clim_files$climSd, paste0(out.dir, '/UVAFME2018_climate_stddev.csv'),
          row.names = FALSE, quote = FALSE)
write.csv(clim_files$extraMN, paste0(out.dir, '/UVAFME2018_climate_ex.csv'),
          row.names = FALSE, quote = FALSE)
write.csv(clim_files$extraSd, paste0(out.dir, 
                                     '/UVAFME2018_climate_ex_stddev.csv'),
          row.names = FALSE, quote = FALSE)
write.csv(lightning_dat, paste0(out.dir, '/UVAFME2018_lightning.csv'), 
          row.names = FALSE, quote = FALSE, na = '')

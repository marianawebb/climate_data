## Data Processing Workflow: Gridmet
# 
# Created by Mari Webb                                                    
# April 27, 2023
# 
# Inputs: gridmet geoTiff, HUC8 watershed shapefiles for CA, NV, UT, AR, CO, NM                                             
# Outputs: CSV of precipitation and temperature trend values for each HUC
# 


# Load Packages -----------------------------------------------------------

library(terra)
library(sf)
library(tidyverse)
library(lubridate)


# Set user-specified variables --------------------------------------------

save_folder <- "data/gridmet/"


# Read in HUC polygons ----------------------------------------------------

# get southwest HUCs
hucfile = paste0('sw_huc8_clean')

# read in HUCX shapefiles
huc_shp = st_read('data/HUCs/', hucfile) %>% 
  mutate(huc8 = as.integer(huc8))


# # get lat long boundaries for SW watersheds, used to crop raster data
# # huc_shp.total_bounds
# huc_ext = ext(st_bbox(huc_shp)['xmin'], st_bbox(huc_shp)['xmax'], 
#             st_bbox(huc_shp)['ymin'], st_bbox(huc_shp)['ymax'])


# Organize climate data netCDFs -------------------------------------------

# get list of netcdf files (in this case, there is just one)
gridmet_list = list.files('data/gridmet/raw/', pattern="", full.names = TRUE)
print(gridmet_list)

gridmet_temp <- terra::rast("data/gridmet/raw/gridmet_trend_temp_F.slope.tif")
gridmet_precip <- terra::rast("data/gridmet/raw/gridmet_trend_precip_in.slope.tif")


# create watersheds mask
huc8_temp = exactextractr::exact_extract(gridmet_temp, 
                                         huc_shp, 
                                         fun = "mean", 
                                         progress = T, 
                                         append_cols = T) %>% 
  dplyr::select(huc8, name, mean) %>% 
  dplyr::rename("temp_trend" = "mean")
huc8_precip = exactextractr::exact_extract(gridmet_precip,
                                           huc_shp, 
                                           fun = "mean", 
                                           progress = T, 
                                           append_cols = T) %>% 
  dplyr::select(huc8, name, mean) %>% 
  dplyr::rename("precip_trend" = "mean")

gridmet_trends <- huc8_temp %>% 
  left_join(huc8_precip, by=c("huc8", "name"))

write_csv(gridmet_trends, "data/gridmet/gridmet_intersect.csv")

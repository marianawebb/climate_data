## Data Post-Processing Workflow: Instances and Trends of Climate Extremes
# 
# Created by Mari Webb                                                    
# March 31, 2023
# 
# Inputs: time series csvs for flood, fire, and drought for each HUC8 in the SW                                            
# Outputs: compounding climate extreme event risk heat map
# # 



# Load Packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(sf)
library(broom)
library(EflowStats)
library(trend)
library(cowplot)
library(stringr)
library(scales)

# Set user-specified variables --------------------------------------------

save_folder <- "data/spei/"


# Read in HUC polygons ----------------------------------------------------

# read in HUC8 shapefiles
huc_shp = st_read('data/HUCs/sw_huc8_clean.shp') %>% 
  mutate(huc8 = as.integer(huc8))

# ggplot() +
#   geom_sf(data = huc_shp, fill="#69b3a2", color="white") +
#   theme_void() 


# read in states shapefile
states_shp <- st_read("data/shapefiles/SW_state_boundaries.shp")


# Organize USDM data netCDFs -------------------------------------------

# get list of netcdf files (in this case, there is just one)
dlist = list.files('data/drought/USDM/', pattern="\\.csv$", full.names = TRUE)
print(dlist)

raw_dta <- lapply(dlist, read.csv) %>% 
  bind_rows()

usdm_df <- raw_dta %>% 
  rowwise() %>% 
  mutate(drought_per = sum(D4,D3,D2,D1),
         drought_12 = sum(D1,D2),
         drought_34 = sum(D3,D4),
         drought_day = ifelse(drought_per> 0.01, 1, 0)) %>% 
  dplyr::rename("huc8"="HUCId", "name"="HUC") %>% 
  mutate(Date = ValidStart,
         WY = get_waterYear(as.Date(Date)))

# ggplot(usdm_df) +
#   geom_histogram(aes(drought_per), bins = 10)
# 
# ggplot(usdm_df) +
#   geom_histogram(aes(drought_day))

usdm_sum <- usdm_df %>% 
  group_by(huc8) %>% 
  summarise(droughts = sum(drought_day))

usdm_mean <- usdm_df %>% 
  group_by(WY, huc8, name) %>% 
  # summarise(annual_drought_per = mean(drought_per, na.rm = T), 
  #           annual_drought_D12 = mean(D1, D2), 
  #           annual_drought_D34 = mean(D3,D4), .groups = "keep") %>% 
  group_by(name, huc8) %>% 
  summarise(mean_drought_per = mean(drought_per, na.rm = T), .groups = "keep")

usdm_per_trends <- usdm_df %>% 
  group_by(WY, huc8, name) %>% 
  # summarise(annual_drought_per = mean(drought_per, na.rm = T), .groups = "keep") %>% 
  group_by(name, huc8) %>% # we group by lon and lat to perform the calculation in each cell
  arrange(WY) %>% 
  summarise(slope = sens.slope(drought_per)$estimates *10,
            sign = mk.test(drought_per)$p.value, .groups = "keep")

usdm_days_trends <- usdm_df %>% 
  group_by(WY, huc8, name) %>% 
  summarise(annual_drought_days = sum(drought_day, na.rm = T), .groups = "keep") %>% 
  group_by(name, huc8) %>% # we group by lon and lat to perform the calculation in each cell
  arrange(WY) %>% 
  summarise(slope_days = sens.slope(annual_drought_days)$estimates *10,
            sign_days = mk.test(annual_drought_days)$p.value, .groups = "keep")

# Visualize Drought -------------------------------------------------------

huc8_drought <- huc_shp %>% 
  left_join(usdm_sum, by = 'huc8') %>% 
  left_join(usdm_mean, by = "huc8") %>% 
  left_join(usdm_per_trends, by = "huc8") %>% 
  left_join(usdm_days_trends, by = "huc8") %>% 
  mutate(slope = ifelse(sign < 0.05, slope, NA),
         slope_days = ifelse(sign_days < 0.05, slope_days, NA))


# map of drought counts in each watershed
ggplot() +
  geom_sf(data = huc8_drought, aes(fill=droughts), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = 1) +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "USDM Drought Days 2000-2022", fill = "# Drought Days")

ggsave("figures/drought/USDM_count_map.png", width = 10, height = 7)

ggplot() +
  geom_sf(data = huc8_drought, aes(fill=slope), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = 1, limits=c(0, 40)) +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "USDM Drought Percent Area Trend 2000-2022", fill = "Sen's Slope")

ggsave("figures/drought/USDM_percent_senslope_map.png", width = 10, height = 7)

ggplot() +
  geom_sf(data = huc8_drought, aes(fill=slope_days), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = 1, limits=c(0, 15)) +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "USDM Annual Drought Days Trend 2000-2022", fill = "Sen's Slope")

ggsave("figures/drought/USDM_days_senslope_map.png", width = 10, height = 7)

ggplot() +
  geom_sf(data = huc8_drought, aes(fill=mean_drought_per), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = 1) +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "USDM Mean Area in Moderate to Exceptional Drought 2000-2022", fill = "% Area in Drought")

ggsave("figures/drought/USDM_mean_percent_map.png", width = 10, height = 7)


# Organize SPEI data netCDFs -------------------------------------------

# # get list of netcdf files (in this case, there is just one)
# dlist = list.files('data/drought/SPEI/spei_huc8s', pattern="\\.csv$", full.names = TRUE)
# print(dlist)
# 
# raw_dta <- lapply(dlist, read.csv) %>% 
#   bind_rows()
# 
# spei_df <- raw_dta %>% 
#   plyr::rename(c('X.4'='D4', 'X.3'='D3', 'X.2'='D2','X.1'='D1', 'X0'='DW0', 'X1'='W1', 'X2'='W2', 'X3'='W3', 'X4'='W4')) %>% 
#   replace(is.na(.), 0) %>% 
#   # mutate(tot_pixels = D4+D3+D2+D1+DW0+W1+W2+W3+W4) %>% 
#   rowwise() %>% 
#   mutate(tot_pixels = sum(D4+D3+D2+D1+DW0+W1+W2+W3+W4),
#          D1_4 = sum(D4+D3+D2+D1),
#          D1_4_per = D1_4/tot_pixels,
#          drought_day = ifelse(D1_4_per> 0.5, 1, 0))
# 
# 
# ggplot(spei_df) +
#   geom_histogram(aes(D1_4_per), bins = 10)
# 
# ggplot(spei_df) +
#   geom_histogram(aes(drought_day))
# 
# spei_sum <- spei_df %>% 
#   group_by(id) %>% 
#   summarise(drought_days = sum(drought_day, na.rm = T))
# 
# spei_events <- spei_df %>% 
#   group_by(id) %>% 
#   summarise(droughts = )
# 
# spei_mean <- spei_df %>% 
#   mutate(WY = get_waterYear(as.Date(Date))) %>% 
#   group_by(WY, id, name) %>% 
#   summarise(annual_drought_per = mean(D1_4_per, na.rm = T), .groups = "keep") %>% 
#   group_by(name, id) %>% 
#   summarise(mean_drought_per = mean(annual_drought_per, na.rm = T))
# 
# spei_trends <- spei_df %>% 
#   mutate(WY = get_waterYear(as.Date(Date))) %>% 
#   group_by(WY, id, name) %>% 
#   summarise(annual_drought_per = mean(D1_4_per, na.rm = T), .groups = "keep") %>% 
#   group_by(name, id) %>% # we group by lon and lat to perform the calculation in each cell
#   summarise(slope = sens.slope(annual_drought_per)$estimates *10,
#             sign = mk.test(annual_drought_per)$p.value, .groups = "keep")
# 
# 
# 
# # Visualize Drought -------------------------------------------------------
# 
# huc8_drought <- huc_shp %>% 
#   left_join(spei_sum, by = c("huc8" = "id")) %>% 
#   left_join(spei_mean, by = c("huc8" = "id")) %>% 
#   left_join(spei_trends, by = c("huc8" = "id")) %>% 
#   mutate(slope = ifelse(sign < 0.05, slope, NA))
# 
# 
# # map of drought counts in each watershed
# ggplot() +
#   geom_sf(data = huc8_drought, aes(fill=drought_days), color="black", lwd = 0.1) +
#   theme_void() +
#   scale_fill_viridis_c(direction = 1) +
#   geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
#   theme_cowplot() +
#   labs(title = "SPEI Drought Count 1980-2022", fill = "# Droughts")
# 
# ggsave("figures/drought/SPEI_count_map.png", width = 10, height = 7)
# 
# ggplot() +
#   geom_sf(data = huc8_drought, aes(fill=slope), color="black", lwd = 0.1) +
#   theme_void() +
#   scale_fill_viridis_c(direction = 1) +
#   scale_color_viridis_c(direction = 1) +
#   geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
#   theme_cowplot() +
#   labs(title = "SPEI Drought Trend 1980-2022", fill = "Sen's Slope")
# 
# ggsave("figures/drought/SPEI_senslope_map.png", width = 10, height = 7)
# 
# ggplot() +
#   geom_sf(data = huc8_drought, aes(fill=mean_drought_per), color="black", lwd = 0.1) +
#   theme_void() +
#   scale_fill_viridis_c(direction = 1) +
#   scale_color_viridis_c(direction = 1) +
#   geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
#   theme_cowplot() +
#   labs(title = "SPEI Mean Area in Moderate to Exceptional Drought 1980-2022", fill = "% Area in Drought")
# 
# ggsave("figures/drought/SPEI_mean_percent_map.png", width = 10, height = 7)



# Fire --------------------------------------------------------------------

# get list of netcdf files (in this case, there is just one)
fire_df = read_csv("data/fire/MTBS/Burn_Intersect.csv") %>% 
  filter(Incid_Type == "Wildfire")

fire_sum <- fire_df %>% 
  group_by(huc8, name) %>% 
  summarise(fire_ct = n())

fire_complete <- fire_df %>% 
  group_by(huc8, name) %>% 
  mutate(Date = as.Date(Ig_Date)) %>% 
  complete(Date = seq.Date(as.Date('1984-01-01'), as.Date('2018-12-31'), by='day')) %>% 
  mutate(burn_percentage = ifelse(is.na(burn_percentage), 0, burn_percentage))

fire_mean <- fire_complete %>% 
  mutate(WY = get_waterYear(Date)) %>% 
  group_by(WY, huc8, name) %>% 
  summarise(annual_burn_per = sum(burn_percentage, na.rm = T), .groups = "keep") %>% 
  group_by(huc8, name) %>% 
  summarise(fire_mean = mean(annual_burn_per))

fire_per_trend <- fire_complete %>% 
  mutate(WY = get_waterYear(Date)) %>% 
  group_by(WY, huc8, name) %>% 
  summarise(annual_burn_per = sum(burn_percentage, na.rm = T), .groups = "keep") %>% 
  group_by(name, huc8) %>% # we group by lon and lat to perform the calculation in each cell
  arrange(WY) %>% 
  summarise(slope = sens.slope(annual_burn_per)$estimates *10,
            sign = mk.test(annual_burn_per)$p.value, .groups = "keep") %>% 
  mutate(slope = ifelse(sign < 0.05, slope, NA))

#no statistically significant increase in the annual count of fires
fire_count_trend <- fire_complete %>% 
  mutate(WY = get_waterYear(Date)) %>% 
  group_by(WY, huc8, name) %>% 
  summarise(annual_burn_ct = n(), .groups = "keep") %>% 
  group_by(name, huc8) %>% # we group by lon and lat to perform the calculation in each cell
  arrange(WY) %>% 
  summarise(slope_ct = sens.slope(annual_burn_ct)$estimates *10,
            sign_ct = mk.test(annual_burn_ct)$p.value, .groups = "keep") %>% 
  mutate(slope_ct = ifelse(sign_ct < 0.05, slope_ct, NA))



# visualize fire instances and trends
huc8_fire <- huc_shp %>% 
  left_join(fire_sum, by = "huc8") %>% 
  left_join(fire_mean, by = "huc8") %>% 
  left_join(fire_per_trend, by = "huc8") %>% 
  mutate(fire_mean = ifelse(is.na(fire_mean), 0, fire_mean),
         fire_ct = ifelse(is.na(fire_ct), 0, fire_ct))

# map of fire counts in each watershed
ggplot() +
  geom_sf(data = huc8_fire, aes(fill=fire_ct), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = -1, option = "F", limits=c(0, 135)) +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "MTBS Fire Events 1984-2021", fill = "# Fires")

ggsave("figures/fire/MTBS_fire_count_map.png", width = 10, height = 7)


# map of mean area burned in each watershed
ggplot() +
  geom_sf(data = huc8_fire, aes(fill=fire_mean), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = -1, option = "F") +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "MTBS Fire Mean Annual Percent Burned 1984-2021", fill = "% Burned")

ggsave("figures/fire/MTBS_fire_mean_map.png", width = 10, height = 7)


# map of burn trends
# map of mean area burned in each watershed
ggplot() +
  geom_sf(data = huc8_fire, aes(fill=slope), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = -1, option = "F") +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "MTBS Fire Annual Trend 1980-2021", fill = "Sen's Slope")

ggsave("figures/fire/MTBS_fire_trend_map.png", width = 10, height = 7)




# Landslide ---------------------------------------------------------------

# get list of netcdf files (in this case, there is just one)
ls_df = read_csv("data/landslides/US_Landslide_1/USGS_landslide_intersect.csv") %>% 
  mutate(yr = str_extract(Date, "\\d{4}")) 

#might want to filter to more recent landslides?
ls_sum <- ls_df %>% 
  filter(yr >= 1900 & yr < 2024) %>% 
  group_by(huc8, name) %>% 
  summarise(ls_ct = n())

ls_mean <- ls_df %>% 
  group_by(huc8, name) %>% 
  summarise(ls_mean = mean(ls_percentage))


# ls_trend_all <- ls_df %>% 
#   dplyr::select(states, huc8, name, Date, lsIntersect_area, ls_percentage) %>% 
#   group_by(huc8, name) %>% 
#   mutate(yr = as.integer(str_extract(Date, "\\d{4}"))) %>% 
#   filter(yr >= 1900 & yr < 2024) %>% 
#   complete(yr = seq(from = 1900, to = 2013, by =1)) %>% 
#   mutate(ls_percentage = ifelse(is.na(ls_percentage), 0, ls_percentage))
# 
# ls_trend <- ls_trend_all %>% 
#   group_by(yr, huc8, name) %>% 
#   summarise(annual_ls_per = sum(ls_percentage, na.rm = T), .groups = "keep") %>% 
#   group_by(name, huc8) %>% # we group by lon and lat to perform the calculation in each cell
#   summarise(slope = sens.slope(annual_ls_per)$estimates *10,
#             sign = mk.test(annual_ls_per)$p.value, .groups = "keep") %>% 
#   mutate(slope = ifelse(sign < 0.05, slope, NA))

# visualize landslide instances and trends
huc8_ls <- huc_shp %>% 
  left_join(ls_sum, by = "huc8") %>% 
  left_join(ls_mean, by = "huc8") %>% 
  mutate(ls_mean = ifelse(is.na(ls_mean), 0, ls_mean),
         ls_ct = ifelse(is.na(ls_ct), 0, ls_ct))

# map of landslide counts in each watershed
ggplot() +
  geom_sf(data = huc8_ls, aes(fill=ls_ct), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = -1, option = "E", limits=c(0, 4000)) +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "USGS Landslide Events 1900-2013", fill = "# Landslides")

ggsave("figures/landslide/USGS_landslide_count_map.png", width = 10, height = 7)


# map of mean area burned in each watershed
ggplot() +
  geom_sf(data = huc8_fire, aes(fill=fire_mean), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = -1, option = "E") +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "USGS Percent Landslide Events 1900-2013", fill = "% Landslides")

ggsave("figures/landslide/USGS_landslide_mean_map.png", width = 10, height = 7)




# Flood -------------------------------------------------------------------


# get list of netcdf files (in this case, there is just one)
flood_df = read_csv("data/flood/modis_dSWE/dSWE_flooding_intersect.csv") %>% 
  mutate(Date = as.Date(paste0(year,"-",month,"-01")))

flood_sum <- flood_df %>%  
  group_by(huc8, name) %>% 
  mutate(mean_flood = mean(mean),
         flood = ifelse(mean > mean_flood, 1, 0)) %>% 
  dplyr::filter(flood == 1) %>% 
  summarise(flood_ct = n())

flood_mean <- flood_df %>% 
  group_by(huc8, name) %>% 
  summarise(flood_mean = mean(mean)*100)

flood_trend <- flood_df %>% 
  mutate(WY = get_waterYear(Date)) %>% 
  group_by(WY, huc8, name) %>% 
  summarise(annual_flood_per = mean(mean, na.rm = T)*100, .groups = "keep") %>% 
  group_by(name, huc8) %>% # we group by lon and lat to perform the calculation in each cell
  arrange(WY) %>% 
  summarise(slope = sens.slope(annual_flood_per)$estimates *10,
            sign = mk.test(annual_flood_per)$p.value, .groups = "keep") %>% 
  mutate(slope = ifelse(sign < 0.05, slope, NA))

# visualize flood instances and trends
huc8_flood <- huc_shp %>% 
  left_join(flood_sum, by = "huc8") %>% 
  left_join(flood_mean, by = "huc8") %>% 
  left_join(flood_trend, by = "huc8") %>% 
  mutate(flood_mean = ifelse(is.na(flood_mean), 0, flood_mean),
         flood_ct = ifelse(is.na(flood_ct), 0, flood_ct))

# map of flood counts in each watershed
ggplot() +
  geom_sf(data = huc8_flood, aes(fill=flood_ct), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = -1, option = "G", limits=c(45, 150)) +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "MODIS Surface Water Extent Flood Events 2003-2019", fill = "# Floods")

ggsave("figures/flood/MODIS_flood_count_map.png", width = 10, height = 7)


# map of mean area flooded in each watershed
ggplot() +
  geom_sf(data = huc8_flood, aes(fill=flood_mean), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = -1, option = "G") +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "MODIS Surface Water Extent Mean Percent 2003-2019", fill = "% Flooded")

ggsave("figures/flood/MODIS_flood_mean_map.png", width = 10, height = 7)


# map of flood trend in each watershed
ggplot() +
  geom_sf(data = huc8_flood, aes(fill=slope), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = -1, option = "G") +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "MODIS Surface Water Extent Annual Trend 2003-2019", fill = "Sen's Slope")

ggsave("figures/flood/MODIS_flood_trend_map.png", width = 10, height = 7)




# # get list of netcdf files (in this case, there is just one)
# flood_df = read_csv("data/flood/DFO/Flood_Intersect.csv") %>% 
#   mutate(Date = as.Date(BEGAN))

# flood_sum <- flood_df %>%  
#   group_by(huc8, name) %>% 
#   summarise(flood_ct = n())
# 
# flood_complete <- flood_df %>% 
#   group_by(huc8, name, Date) %>%
#   rowwise() %>%
#   do(data.frame(huc8=.$huc8, name=.$name, start=.$BEGAN, end=.$ENDED, cause=.$MAINCAUSE, 
#                 severity=.$SEVERITY, floodIntersect_area=.$floodIntersect_area,
#                 flood_percentage=.$flood_percentage, Date=seq(.$BEGAN,.$ENDED,by="1 day"))) %>% 
#   group_by(huc8, name) %>% 
#   complete(Date = seq.Date(as.Date('1985-01-01'), as.Date('2021-12-31'), by='day')) %>% 
#   mutate(flood_percentage = ifelse(is.na(flood_percentage), 0, flood_percentage))

# flood_mean <- flood_complete %>% 
#   group_by(huc8, name) %>% 
#   summarise(flood_mean = mean(flood_percentage))
# 
# flood_trend <- flood_complete %>% 
#   mutate(WY = get_waterYear(Date)) %>% 
#   group_by(WY, huc8, name) %>% 
#   summarise(annual_flood_per = mean(flood_percentage, na.rm = T), .groups = "keep") %>% 
#   group_by(name, huc8) %>% # we group by lon and lat to perform the calculation in each cell
#   arrange(WY) %>% 
#   summarise(slope = sens.slope(annual_flood_per)$estimates *10,
#             sign = mk.test(annual_flood_per)$p.value, .groups = "keep") %>% 
#   mutate(slope = ifelse(sign < 0.05, slope, NA))


# # visualize flood instances and trends
# huc8_flood <- huc_shp %>% 
#   left_join(flood_sum, by = "huc8") %>% 
#   left_join(flood_mean, by = "huc8") %>% 
#   mutate(flood_mean = ifelse(is.na(flood_mean), 0, flood_mean),
#          flood_ct = ifelse(is.na(flood_ct), 0, flood_ct))
# 
# # map of flood counts in each watershed
# ggplot() +
#   geom_sf(data = huc8_flood, aes(fill=flood_ct), color="black", lwd = 0.1) +
#   theme_void() +
#   scale_fill_viridis_c(direction = -1, option = "G", limits=c(0, 20)) +
#   geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
#   theme_cowplot() +
#   labs(title = "DFO Flood Events 1985-2021", fill = "# Floods")
# 
# ggsave("figures/flood/DFO_flood_count_map.png", width = 10, height = 7)
# 
# 
# # map of mean area flooded in each watershed
# ggplot() +
#   geom_sf(data = huc8_flood, aes(fill=flood_mean), color="black", lwd = 0.1) +
#   theme_void() +
#   scale_fill_viridis_c(direction = -1, option = "G",) +
#   geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
#   theme_cowplot() +
#   labs(title = "DFO Mean Percent Flooded 1985-2021", fill = "% Flooded")
# 
# ggsave("figures/flood/DFO_flood_mean_map.png", width = 10, height = 7)





# Create normalized and scaled heatmap ------------------------------------

huc8_normalized <- flood_mean %>% 
  left_join(usdm_mean, by = c("huc8", "name")) %>% 
  left_join(fire_mean, by = c("huc8", "name")) %>% 
  left_join(flood_trend %>% 
              dplyr::rename("flood_slope"="slope") %>% 
              dplyr::filter(sign<0.05), 
            by = c("huc8", "name")) %>% 
  left_join(usdm_per_trends %>% 
              dplyr::rename("usdm_slope"="slope") %>% 
              dplyr::filter(sign<0.05),
            by = c("huc8", "name")) %>% 
  left_join(fire_per_trend %>% 
              dplyr::rename("fire_slope"="slope") %>% 
              dplyr::filter(sign<0.05),
            by = c("huc8", "name")) %>% 
  ungroup() %>% 
  mutate_at(c(3:5), ~replace_na(.,0)) %>% ## NEED TO CHECK ON SOME NA VALUES INCLUDING AGUA FRIA WATERSHED, also need to figure out what to do with NA values, they should probably be assigned 0 for fire
  # mutate(across(3:5, ~ (c(scale(., scale=F))))) %>% 
  mutate_at(c(3:5), ~rescale(c(scale(.)), to = c(0,1))) %>%
  dplyr::rename("drought_mean" = "mean_drought_per") %>% 
  mutate(flood_ct = ifelse(flood_slope > 0, 1, 0),
         usdm_ct = ifelse(usdm_slope > 0, 1, 0),
         fire_ct = ifelse(fire_slope > 0, 1, 0)) %>% 
  rowwise() %>% 
  mutate(risk_index = (mean(c(drought_mean, flood_mean, fire_mean))),
         name = ifelse(name == "Aqua Fria", "Agua Fria", name),
         risk_ct = sum(flood_ct, usdm_ct, fire_ct, na.rm = T))
         

huc8_risk <- huc_shp %>% 
  left_join(huc8_normalized, by = c("huc8", "name")) %>% 
  mutate(risk_index = rescale(risk_index, to= c(0,5)))

huc8_hist <- huc8_normalized %>% pivot_longer(cols = 3:5, names_to = "names")

ggplot(huc8_normalized %>% pivot_longer(cols = 3:5, names_to = "names"), aes(value)) +
  geom_histogram() +
  facet_wrap(~names)


# map of risk in each watershed
ggplot() +
  geom_sf(data = huc8_risk, aes(fill=risk_index), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_viridis_c(direction = 1, option = "H", limits = c(0,5)) +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "Climate Extreme Events Heat Map", fill = "Compounding\nRisk Index")

ggsave("figures/heatmap/risk_index_map.png", width = 10, height = 7)


ggplot() +
  geom_sf(data = huc8_risk, aes(fill=risk_ct), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_distiller(palette = "Purples", direction = 1) +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "Climate Extreme Events Trend Map", fill = "# Increasing\nExtremes")

 ggsave("figures/heatmap/compounding_trends_map.png", width = 10, height = 7)


# Gridmet Trends Map ------------------------------------------------------

gridmet_trends <- read_csv("data/gridmet/gridmet_intersect.csv")

# visualize flood instances and trends
huc8_gridmet <- huc_shp %>% 
  left_join(gridmet_trends, by = "huc8")



# map of flood trend in each watershed
ggplot() +
  geom_sf(data = huc8_gridmet, aes(fill=temp_trend), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_distiller(direction = 1, palette = "OrRd") +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "Gridmet Temperature (F) Trend 1980-2022", fill = "Sen's Slope")

ggsave("figures/gridmet/gridmet_temp_trends_map.png", width = 10, height = 7)


# map of flood trend in each watershed
ggplot() +
  geom_sf(data = huc8_gridmet, aes(fill=precip_trend), color="black", lwd = 0.1) +
  theme_void() +
  scale_fill_distiller(direction = 1, palette = "GnBu", limits = c(-0.5,0.5)) +
  # scale_fill_viridis_c(direction = -1, option = "E", 
  #                      na.value = "lightgray") +
  geom_sf(data = states_shp, fill = NA, lwd = 0.5, color = "black") +
  theme_cowplot() +
  labs(title = "Gridmet Precip (in.) Trend 1980-2022", fill = "Sen's Slope") 

ggsave("figures/gridmet/gridmet_precip_trends_map.png", width = 10, height = 7)


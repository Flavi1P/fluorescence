library(tidyverse)
library(ncdf4)
library(readxl)
library(sf)
library(oce)
library(metR)

# soclim <- read_delim("Data/soclim/soclim_btl_V9.txt", 
#            "\t", escape_double = FALSE, trim_ws = TRUE)
# soclim <- janitor::clean_names(soclim[, c(1:103)])


soclim <- read_csv("Data/soclim_abs.csv")

soclim <- soclim %>% mutate(cast_id = paste(station, ctd_number))

soclim$total_chlorophyll_a <- as.numeric(soclim$total_chlorophyll_a)
ggplot(soclim)+
  geom_point(aes(x = fluorescence_eco_rfu, y = - depth_m, group = cast_id))+
  geom_point(aes(x = fluorescence_wetlabs_rfu, y = - depth_m), colour = "Brown")+
  geom_point(aes(x = total_chlorophyll_a, y = -depth_m), colour = "Green")+
  ylim(-400,0)+
  facet_wrap(.~cast_id)

ggplot(soclim)+
  geom_point(aes(x = longitude_deg_east, y = latitude_deg_north, colour = station))

soclim_ratio <- soclim %>% na.omit() %>%
  mutate(relativ_chl = total_chlorophyll_a/max(total_chlorophyll_a),
                                  ratio = fluorescence_eco_rfu*2.5/total_chlorophyll_a,
                                  relativ_ratio = fluorescence_eco_rfu/relativ_chl,
         yield = fluorescence_eco_rfu/pur470) %>% 
  filter(depth_m < 100) %>% 
  group_by(cast_id) %>% 
  select(cast_id, "lon" = longitude_deg_east, "lat" = latitude_deg_north, relativ_chl, ratio, relativ_ratio, total_chlorophyll_a, fluorescence_eco_rfu, yield) %>% 
  summarise_all(mean, na.omit = TRUE)

worldmap <- map_data("world") %>% filter(long > min(soclim_ratio$lon) & long < max(soclim_ratio$lon) & lat < max(soclim_ratio$lat) & lat > min(soclim_ratio$lat))

ggplot(soclim_ratio)+
  geom_point(aes(x = lon, y = lat, colour = ratio), size = 4)+
  geom_polygon(aes(x = long, y = lat), data = worldmap)+
  scale_colour_distiller(palette = "RdBu")+
  coord_quickmap()+
  theme_minimal()

soclim_plot <- soclim %>% na.omit() %>%
  mutate(relativ_chl = total_chlorophyll_a/max(total_chlorophyll_a),
         ratio = fluorescence_eco_rfu*2.5/total_chlorophyll_a,
         relativ_ratio = fluorescence_eco_rfu/relativ_chl) %>% 
  filter(depth_m < 100) %>% 
  mutate(group = case_when(
    ratio < 1.9 ~ "low",
    ratio >= 1.9 ~ "high"))

ggplot(soclim_plot)+
  geom_point(aes(x = total_chlorophyll_a, y = fluorescence_eco_rfu, colour = group))+
  theme_bw()

ggplot(soclim_plot)+
  geom_point(aes(x = longitude_deg_east, y = latitude_deg_north, colour = group), size = 4)+
  geom_polygon(aes(x = long, y = lat), data = worldmap)+
  coord_quickmap()+
  theme_minimal()

ggplot(soclim_plot)+
  geom_point(aes(y = -depth_m, x = latitude_deg_north, colour = group), size = 4)

soclim_dist <- soclim %>% filter(fluorescence_eco_rfu > 0) %>%
  mutate(log_chla = log(total_chlorophyll_a),
                            log_eco = log(fluorescence_eco_rfu),
         yield = fluorescence_eco_rfu/pur470)
ggplot(soclim_dist)+
  geom_density(aes(x = log_chla), colour = "Brown")+
  geom_density(aes(x = log_eco), colour = "Grey")+
  theme_bw()

ggplot(soclim_dist)+
  geom_point(aes(x = fluorescence_eco_rfu, y = - depth_m, group = cast_id))+
  geom_point(aes(x = fluorescence_eco_rfu, y = -depth_m), data = filter(soclim_dist, log_eco < -4), colour = "Red", size = 2)+
  geom_point(aes(x = total_chlorophyll_a, y = -depth_m), colour = "Green")+
  ylim(-400,0)+
  facet_wrap(.~cast_id)

soclim_dist_clean <- soclim_dist %>% filter(log_eco > -4)

ggplot(soclim_dist_clean)+
  geom_density(aes(x = log_chla), colour = "Brown")+
  geom_density(aes(x = log_eco), colour = "Grey")+
  theme_bw()

ggplot(soclim_dist_clean)+
  geom_point(aes(x =log_chla, y = log_eco, colour = station))

ggplot(soclim_dist_clean)+
  geom_point(aes(x = total_chlorophyll_a, y = fluorescence_eco_rfu, colour = station))

soclim_dist_clean <- soclim_dist_clean %>%
  na.omit() %>%
  mutate(relativ_chl = total_chlorophyll_a/max(total_chlorophyll_a),
         ratio = fluorescence_eco_rfu*2.5/total_chlorophyll_a,
         relativ_ratio = fluorescence_eco_rfu/relativ_chl,
         eco_test = fluorescence_eco_rfu * max(total_chlorophyll_a)) %>%
  mutate(group = case_when(
    ratio < 1.9 ~ "low",
    ratio >= 1.9 ~ "high"))

ggplot(soclim_dist_clean)+
  geom_point(aes(x = total_chlorophyll_a, y = fluorescence_eco_rfu, colour = group))

ggplot(filter(soclim_dist_clean, group == "high"))+
  geom_point(aes(x = eco_test, y = - depth_m, group = cast_id))+
  geom_point(aes(x = total_chlorophyll_a, y = -depth_m), colour = "Green")+
  ylim(-200,0)+
  facet_wrap(.~cast_id)

ggplot(soclim_dist_clean)+
  geom_point(aes(x = eco_test, y = - depth_m, colour = group))+
  geom_point(aes(x = total_chlorophyll_a, y = -depth_m), colour = "Green")+
  ylim(-200,0)+
  facet_wrap(.~cast_id, scales = 'free_x')

ggplot(soclim_dist_clean)+
  geom_point(aes(x = yield, y =  -depth_m, colour = group))+
  ylim(-200,0)+
  facet_wrap(.~cast_id)

ggplot(soclim_dist_clean)+
  geom_violin(aes(x = group, y = yield))+
  geom_jitter(aes(x = group, y = yield))+
  theme_bw()
  
ggplot(soclim_dist_clean)+
  geom_violin(aes(x = group, y = ratio))+
  geom_jitter(aes(x = group, y = ratio))+
  theme_bw()

nc_path <- c("Data/soclim/current.nc")

nc <- nc_open(nc_path)

north_mat <- ncvar_get(nc, "vo")
east_mat <- ncvar_get(nc, "uo")
longitude <- ncvar_get(nc, "longitude")
latitude <- ncvar_get(nc, "latitude")

north_df <- as.data.frame(north_mat)
names(north_df) <- latitude

east_df <- as.data.frame(east_mat)
names(east_df) <- latitude

north_df <- north_df %>% mutate("lon"=  longitude) %>% 
  pivot_longer(c(1:601), names_to = "lat", values_to = "north")

east_df <- east_df %>% mutate("lon" = longitude) %>% 
  pivot_longer(c(1:601), names_to = "lat", values_to = "east")

current_df <- left_join(north_df, east_df) %>%
  mutate(vel = sqrt(north^2+east^2),
         lat = as.numeric(lat))

subset_df <- current_df %>% filter(lon > 50 & lon < 80 & lat < -25 & lat > -60)

mean_df <- subset_df %>% mutate(lon_ron = round(lon),
                                lat_ron = round(lat)) %>% 
  group_by(lon_ron, lat_ron) %>% 
  summarise_all(mean)

# ggplot(mean_df)+
#   geom_raster(data = subset_df, aes(x = lon, y = lat, fill = vel), interpolate = TRUE)+
#   geom_segment(aes(x = lon, xend = lon + east/1.2, y = lat, yend = lat+north/1.2), 
#                arrow = arrow(angle = 20, length = unit(.2, "cm"), type = "open"))+
#   geom_polygon(aes(x = long, y = lat), data = worldmap, fill = "Brown")+
#   scale_fill_gradientn(name = "Speed\n(m/s)",colours = oce::oceColorsVelocity(120), 
  #                      limits = c(0,1.6), breaks =seq(0.1,1.6,.3))+
  # geom_point(aes(x = longitude_deg_east, y = latitude_deg_north, colour = group), data =soclim_dist_clean)+
  # theme_bw()+
  # theme(legend.position = "right",
  #       legend.key.height = unit(1.4, "cm"), 
  #       legend.background = element_blank(),
  #       axis.text = element_text(size = 12, colour = 1))+
  # labs(x = "", y = "")+
  # ggtitle("Current velocity near KI in October 2016")+
  # coord_quickmap()

ggplot() +
  metR::geom_contour_fill(data = subset_df, aes(x = lon, y = lat, z = vel), na.fill = TRUE, bins = 70) + 
  metR::geom_vector(data = mean_df, aes(x = lon, y = lat, dx = east, dy = north), 
                    arrow.angle = 30, arrow.type = "open", arrow.length = .5, 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw")+
  geom_polygon(aes(x = long, y = lat), data = worldmap, fill = "Grey")+
  geom_point(aes(x = longitude_deg_east, y = latitude_deg_north, colour = group), data =soclim_dist_clean)+
  scale_fill_gradientn(name = "Speed\n(m/s)",colours = oceColorsVelocity(120), 
                       limits = c(0,1.5), breaks =seq(0.1,1.6,.3))+
  theme_bw()+
  theme(legend.position = "right",
        legend.key.height = unit(1.2, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1))+
  scale_mag(max = 1, name = "Speed", max_size = 0.5)+
  scale_colour_brewer(palette = "Dark2")+
  labs(x = "", y = "")+
  coord_quickmap()








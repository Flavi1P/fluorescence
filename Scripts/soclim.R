library(tidyverse)
library(ncdf4)
library(readxl)
library(sf)
library(oce)
library(metR)
library(vegan)

# soclim <- read_delim("Data/soclim/soclim_btl_V9.txt", 
#            "\t", escape_double = FALSE, trim_ws = TRUE)
# soclim <- janitor::clean_names(soclim[, c(1:103)])


soclim <- read_csv("Data/soclim/soclim_abs.csv")

soclim <- soclim %>% mutate(cast_id = paste(station, ctd_number))

soclim$total_chlorophyll_a <- as.numeric(soclim$total_chlorophyll_a)
ggplot(soclim)+
  geom_point(aes(x = fluorescence_eco_rfu, y = - depth_m, group = cast_id))+
  geom_point(aes(x = fluorescence_wetlabs_rfu, y = - depth_m), colour = "Brown")+
  geom_point(aes(x = total_chlorophyll_a, y = -depth_m), colour = "Green")+
  ylim(-400,0)+
  facet_wrap(.~cast_id)

soclim_plot <- filter(soclim, fluorescence_eco_rfu != "NA")
ggplot(filter(soclim_plot, !cast_id %in% c("BDT1 7", "O11 11", "O22 1")))+
  geom_point(aes(x = fluorescence_eco_rfu, y = - depth_m, colour = "Eco_RFU"))+
  geom_point(aes(x = total_chlorophyll_a, y = -depth_m, colour = "Chla HPLC"))+
  ylim(-400,0)+
  theme_bw()+
  facet_wrap(.~cast_id)+
  scale_colour_brewer(palette = "Dark2")

#ggsave("Output/soclim/cast_raw.png")

soclim_ratio <- soclim %>% na.omit() %>%
  mutate(relativ_chl = total_chlorophyll_a/max(total_chlorophyll_a),
                                  ratio = fluorescence_eco_rfu/total_chlorophyll_a,
                                  relativ_ratio = fluorescence_eco_rfu/relativ_chl,
         yield = fluorescence_eco_rfu/pur470) %>% 
  filter(depth_m < 100) %>% 
  group_by(cast_id) %>% 
  dplyr::select(cast_id, "lon" = longitude_deg_east, "lat" = latitude_deg_north, relativ_chl, ratio, relativ_ratio, total_chlorophyll_a, fluorescence_eco_rfu, yield) %>% 
  summarise_all(mean, na.omit = TRUE)

worldmap <- map_data("world") %>% filter(long > min(soclim_ratio$lon) & long < max(soclim_ratio$lon) & lat < max(soclim_ratio$lat) & lat > min(soclim_ratio$lat))

ggplot(soclim_ratio)+
  geom_point(aes(x = lon, y = lat, colour = ratio), size = 4)+
  geom_polygon(aes(x = long, y = lat), data = worldmap)+
  scale_colour_distiller(palette = "RdBu", name = "Eco RFU / Tchla HPLC")+
  coord_quickmap()+
  theme_minimal()

#ggsave("Output/soclim/ratio_map.png")

soclim_plot <- soclim %>% filter(fluorescence_eco_rfu != "NA" & total_chlorophyll_a != "NA") %>% 
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
  geom_density(aes(x = log_chla, colour = "Chla HPLC"))+
  geom_density(aes(x = log_eco, colour = "Eco RFU"))+
  theme_bw()

#ggsave("Output/soclim/distribution_double.png")

ggplot(soclim_dist)+
  geom_point(aes(x = fluorescence_eco_rfu, y = - depth_m, colour = "Eco RFU"))+
  geom_point(aes(x = total_chlorophyll_a, y = -depth_m, colour = "Chla HPLC"))+
  geom_point(aes(x = fluorescence_eco_rfu, y = -depth_m), data = filter(soclim_dist, log_eco < -4), colour = "Red", size = 2)+
  ylim(-400,0)+
  facet_wrap(.~cast_id)+
  theme_bw()+
  scale_colour_brewer(palette = "Dark2")

#ggsave(("Output/soclim/cut_cast.png"))

soclim_dist_clean <- soclim_dist %>% filter(log_eco > -4)

ggplot(soclim_dist_clean)+
  geom_density(aes(x = log_chla, colour ="Eco RFU"))+
  geom_density(aes(x = log_eco, colour = "Chla HPLC"))+
  theme_bw()

#ggsave("Output/soclim/distrib_norm.png")

ggplot(soclim_dist_clean)+
  geom_point(aes(x =log_chla, y = log_eco, colour = station))

ggplot(soclim_dist_clean)+
  geom_point(aes(x = total_chlorophyll_a, y = fluorescence_eco_rfu, colour = station))+
  theme_bw()

#ggsave("Output/soclim/raw_scatterplot.png")
# ncdf bathy --------------------------------------------------------------

nc_bath <- "Data/GEBCO_2014_6x6min_Global.nc"

bathy <- nc_open(nc_bath)
lon <- ncvar_get(bathy, "lon")
lat <- ncvar_get(bathy, "lat")
height <- ncvar_get(bathy, "Height")

bath_df <- as.data.frame(height)
names(bath_df) <- lat

bath_df <- bath_df %>% mutate("lon" = lon) %>%
  pivot_longer(c(1:1800), names_to = "lat", values_to = "height") %>% 
  mutate(lat = as.numeric(lat)) %>% 
  filter(lon > 50 & lon < 80 & lat < -25 & lat > -60)
# plot of bathy -----------------------------------------------------------
soclim_station <- dplyr::select(soclim, longitude_deg_east, latitude_deg_north, station) %>%
  group_by(station) %>% 
  summarise_all(mean) %>% 
  ungroup()

ggplot() +
  geom_tile(aes(x = lon, y = lat, fill = height), data = bath_df)+
  geom_polygon(aes(x = long, y = lat), data = worldmap, fill = "Grey")+
  geom_text(aes(x = longitude_deg_east, y = latitude_deg_north, label = station), data =soclim_station)+
  scale_fill_gradientn(name = "Bathymetri",colours = oceColorsGebco(),
                       limits = c(-6000,0), breaks =seq(-5000, 0, 1000))+
  scale_colour_brewer(name = "fluorescence yield", palette = "Set1")+
  theme_bw()+
  theme(legend.position = "right",
        legend.key.height = unit(1.2, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1))+
  labs(x = "", y = "")+
  ggtitle("Soclim stations")+
  coord_quickmap()

#ggsave("Output/soclim/soclim_station.png")

ggplot() +
  geom_tile(aes(x = lon, y = lat, fill = height), data = bath_df)+
  geom_polygon(aes(x = long, y = lat), data = worldmap, fill = "Grey")+
  metR::geom_vector(data = mean_df, aes(x = lon, y = lat, dx = east, dy = north), 
                    arrow.angle = 30, arrow.type = "open", arrow.length = .5, 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw")+
  geom_point(aes(x = longitude_deg_east, y = latitude_deg_north, colour = group), data =soclim_dist_clean)+
  scale_fill_gradientn(name = "Bathymetri",colours = oceColorsGebco(),
                       limits = c(-6000,0), breaks =seq(-5000, 0, 1000))+
  scale_colour_brewer(name = "fluorescence yield", palette = "Set1")+
  scale_mag(max = 1, name = "Speed", max_size = 0.5)+
  theme_bw()+
  theme(legend.position = "right",
        legend.key.height = unit(1.2, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1))+
  labs(x = "", y = "")+
  ggtitle("Soclim stations")+
  coord_quickmap()





# ncdf current ------------------------------------------------------------

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

# resume data -------------------------------------------------------------



mean_df <- subset_df %>% mutate(lon_ron = round(lon),
                                lat_ron = round(lat)) %>% 
  group_by(lon_ron, lat_ron) %>% 
  summarise_all(mean)

# heavy plot of current ---------------------------------------------------


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













# pigment analysis --------------------------------------------------------

soclim_pig <- soclim_dist_clean %>% dplyr::select(station, 'lat' = latitude_deg_north,
                                       'lon' = longitude_deg_east,
                                      'depth' = depth_m,
                                      cast_id,
                                      'fuco' = fucoxanthin,
                                      'zea' = zeaxanthin,
                                      'allo' = alloxanthin,
                                      'peri' = peridinin,
                                      'but' = x19_butanoyloxyfucoxanthin,
                                      'hex' = x19_hexanoyloxyfucoxanthin,
                                      'tchlb' = t_chlb,
                                      'tchla' = total_chlorophyll_a)

soclim_pig[soclim_pig == "LOD"] <- "0"
  
soclim_pig <- soclim_pig %>% na.omit() %>% 
  mutate(fuco = as.numeric(fuco),
         zea = as.numeric(zea),
         allo = as.numeric(allo),
         peri = as.numeric(peri),
         but = as.numeric(but),
         hex = as.numeric(hex),
         tchlb = as.numeric(tchlb),
         wdp = 1.56 * fuco + 0.92 * peri + 4.8 * allo + 1.02 * but + 1.12 * hex + 1.51 * zea + 0.69 * tchlb,
         micro = (1.56 * fuco + 0.92 * peri)/wdp,
         nano = (4.8 * allo + 1.02 * but + 1.51 * hex)/wdp,
         pico = (1.51 * zea + 0.69 * tchlb)/wdp)

#♠perform CA
AFC <- cca(dplyr::select(soclim_pig, fuco:tchlb), scale = TRUE)

#extract scores of eache point
scores <- data.frame(scores(AFC, choices = c(1,2,3), display = "site"))
soclim_pig <- bind_cols(soclim_pig, scores) #addthem to lov df

#compute score of environmental variables a.k.a pigments
fitscore <- envfit(AFC, dplyr::select(soclim_pig, fuco:tchlb), na.rm = TRUE)  
fitarrow <- as.data.frame(fitscore$vectors$arrows)

#plot the all
ggplot(soclim_pig)+
  geom_point(aes(x = CA1, y = CA2, colour = ))+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = fitarrow)+
  geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(fitarrow), fontface = 2), data = fitarrow)+
  scale_color_viridis_c()

#create tree cluster
distsoclim<- dist(dplyr::select(soclim_pig, CA1, CA2))
plot(hclust(distsoclim))
ggsave("Output/soclim/dendrogram.png")

soclim_pig$group <- as.factor(cutree(hclust(distsoclim, method = "ward.D"), k = 3))

ggplot(soclim_pig)+
  geom_point(aes(x = CA1, y = CA2, colour = group))+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = fitarrow)+
  geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(fitarrow), fontface = 2), data = fitarrow)+
  scale_color_brewer(palette = "Set1")

#ggsave("Output/soclim/afc_clust.png")

ggplot(soclim_pig)+
  geom_point(aes(x = lat, y = -depth, colour = group))+
  scale_color_brewer(palette = "Set1")

ggplot() +
  geom_tile(aes(x = lon, y = lat, fill = height), data = bath_df)+
  geom_polygon(aes(x = long, y = lat), data = worldmap, fill = "Grey")+
  geom_point(aes(x = lon, y = lat, colour = group), data =soclim_pig)+
  scale_fill_gradientn(name = "Bathymetri",colours = oceColorsGebco(),
                       limits = c(-6000,0), breaks =seq(-5000, 0, 1000))+
  scale_colour_brewer(name = "fluorescence yield", palette = "Set1")+
  theme_bw()+
  theme(legend.position = "right",
        legend.key.height = unit(1.2, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1))+
  labs(x = "", y = "")+
  ggtitle("Soclim stations")+
  coord_quickmap()

ggsave("Output/soclim/cluster_location.png")

# treemap of pigments -----------------------------------------------------

soclim_pop <- soclim_pig %>% dplyr::select(fuco:tchlb, group) %>% 
  group_by(group) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  pivot_longer(fuco:tchlb, names_to = "pigment", values_to = "concentration")

pig_index <- data.frame("pigment" = c("fuco", "zea", "allo", "peri", "but", "hex", "tchlb"),
                        "size_class" = c("micro", "pico", "nano", "micro", "nano", "nano", "pico"))

soclim_pop <- left_join(soclim_pop, pig_index)
library(treemap)

treemap(filter(soclim_pop, group == "1"), index=c("size_class","pigment"), vSize="concentration",
        
        type="index",                            # How you color the treemap. type help(treemap) for more info
        palette = "Set1",                        # Select your color palette from the RColorBrewer presets or make your own.
        title="Group 1",                      # Customize your title
        fontsize.title=12,                       # Size of the title
        
) 


treemap(filter(soclim_pop, group == "2"), index=c("size_class","pigment"), vSize="concentration",
        
        type="index",                            # How you color the treemap. type help(treemap) for more info
        palette = "Set1",                        # Select your color palette from the RColorBrewer presets or make your own.
        title="Group 2",                      # Customize your title
        fontsize.title=12,                       # Size of the title
        
) 

treemap(filter(soclim_pop, group == "3"), index=c("size_class","pigment"), vSize="concentration",
        
        type="index",                            # How you color the treemap. type help(treemap) for more info
        palette = "Set1",                        # Select your color palette from the RColorBrewer presets or make your own.
        title="Group 3",                      # Customize your title
        fontsize.title=12,                       # Size of the title
        
) 


# yield analysis ----------------------------------------------------------

soclim_clust <- dplyr::select(soclim_pig, group, station, depth)

soclim_dist_clean <- left_join(soclim_dist_clean, soclim_clust, by = c("station", "depth_m" = "depth"))

soclim_dist_clean <- soclim_dist_clean %>%
  na.omit() %>%
  mutate(relativ_chl = total_chlorophyll_a/max(total_chlorophyll_a),
         ratio = fluorescence_eco_rfu/total_chlorophyll_a,
         relativ_ratio = fluorescence_eco_rfu/relativ_chl,
         eco_test = fluorescence_eco_rfu * max(total_chlorophyll_a)) %>%
  mutate(yield_case = case_when(
    ratio < 0.7~ "low",
    ratio >= 0.7 ~ "high")) %>% 
  filter(group != "2")

ggplot(soclim_dist_clean)+
  geom_point(aes(x = total_chlorophyll_a, y = fluorescence_eco_rfu, colour = group))


ggplot(soclim_dist_clean)+
  geom_violin(aes(x = group, y = yield, colour = group))+
  geom_jitter(aes(x = group, y = yield, colour = group))+
  theme_bw()
ggsave("Output/soclim/yield_cluster.png")

ggplot(soclim_dist_clean)+
  geom_violin(aes(x = group, y = ratio, colour = group))+
  geom_jitter(aes(x = group, y = ratio, colour = group))+
  theme_bw()+
  scale_color_brewer(palette = "Set1")

ggplot(soclim_dist_clean)+
  geom_point(aes(x = total_chlorophyll_a, y = fluorescence_eco_rfu, colour = group))+
  geom_smooth(aes(x = total_chlorophyll_a, y = fluorescence_eco_rfu, colour = group), method = "lm", se = FALSE)+
  theme_bw()

model <- aov(yield~group, data = soclim_dist_clean)
summary(model)

fitted_models <-  soclim_dist_clean %>% group_by(group) %>% do(tidy(lm(fluorescence_eco_rfu ~ total_chlorophyll_a, data = .))) %>% ungroup() %>% 
  filter(term != '(Intercept)') %>%  
  select(group, 'estimate_eco' = estimate)

fitted_models <-  soclim_dist_clean %>% group_by(group) %>% do(tidy(lm(eco_test ~ total_chlorophyll_a, data = .))) %>% ungroup() %>% 
  filter(term != '(Intercept)') %>%  
  select(group, 'estimate_eco' = estimate)
plot(TukeyHSD(model))


library(tidyverse)
library(sf)
library(geosphere)
library(tidybgc)

argo <- read_csv("Data/argo_db.csv")
hplc <- read_csv("Data/hplc_argo_campaign")
map <- map_data('world')


#hplc lon lat correction
hplc <- filter(hplc, !is.na(lon) & !is.na(lat))

hplc$lat[hplc$lon>100] <- - hplc$lat[hplc$lon>100]
hplc$lon[hplc$lon>100] <- - hplc$lon[hplc$lon>100]


ggplot()+
  geom_polygon(aes(x = long, y = lat, group= group), data = map)+
  geom_point(aes(x = lon, y = lat), data = argo, size = 4, colour = 'brown')+
  geom_point(aes(x = lon, y = lat), data = hplc, colour = 'grey', alpha = 0.1)+
  coord_quickmap()


hplc <- hplc %>% mutate(profil_id = paste(date, lon, lat, station_name, ctd_number, sep = "_"))
argo <- argo %>% mutate(profil_id = paste(float, lon, lat, sep = "_"))

hplc_lonlat <- select(hplc, lon, lat) %>% distinct() %>% na.omit()
argo_lonlat <- select(argo, lon, lat) %>% distinct()

names(hplc_lonlat) <- c("longitude", "latitude")
names(argo_lonlat) <- c("longitude", "latitude")

match_final <- data.frame("lon_argo" = numeric(),
                          "lat_argo" = numeric(),
                          "lon_hplc" = numeric(),
                          "lat_hplc" = numeric())
for(i in c(1:length(argo_lonlat$longitude))){
  dist <- distGeo(argo_lonlat[i,], hplc_lonlat)
  dist_min <- min(dist)
  hplc_position <- hplc_lonlat[which(dist == dist_min),]
  match <- bind_cols(argo_lonlat[i,], hplc_position)
  names(match) <- c("lon_argo", "lat_argo", "lon_hplc", "lat_hplc")
  match_final <- bind_rows(match_final, match)
}

range(hplc_lonlat$longitude)


ggplot(match_final)+
  geom_polygon(aes(x = long, y = lat, group= group), data = map)+
  geom_point(aes(x = lon_argo, y = lat_argo), size = 4, colour = 'brown')+
  geom_point(aes(x = lon_hplc, y = lat_hplc), colour = 'grey')+
  coord_quickmap()


argo_with_hplc_coord <- argo %>% left_join(match_final, by = c("lon" = "lon_argo", "lat" = "lat_argo"))
argo_with_hplc_data <- argo_with_hplc_coord %>% left_join(hplc, by = c("lon_hplc" = "lon", "lat_hplc" = "lat"))

names(argo_with_hplc_data) <- gsub('.y', '_hplc', colnames(argo_with_hplc_data))
names(argo_with_hplc_data) <- gsub('.x', '_argo', colnames(argo_with_hplc_data))

merged <- argo_with_hplc_data %>% mutate(diff_depth = abs(pres - depth_hplc),
                                                      diff_date = abs(date_argo - date_hplc)) 



ggplot(argo_with_hplc_data)+
  geom_polygon(aes(x = long, y = lat, group= group), data = map)+
  geom_point(aes(x = lon, y = lat), size = 4, colour = 'brown')+
  geom_point(aes(x = lon_hplc, y = lat_hplc), colour = 'grey')+
  coord_quickmap()

for(i in unique(argo$float)){
  data <- filter(argo_with_hplc_data, float == i)
  ggplot(data)+
    geom_point(aes(x = chla_adjusted, y = -pres))+
    geom_point(aes(x = t_chla, y = -depth_hplc), colour = "green")
}

test <- extract_sd(ref_name$path[ref_name$number == 6902737], c("CHLA_ADJUSTED", "PRES"))
ggplot(test)+
  geom_point(aes(y = - pres, x = chla_adjusted))
# Creating new PDF file
pdf("Export_Plots.pdf", width = 16 , height = 10, title = "Argo/HPLC matchup")

# Deciding rows and cols
par(mfrow=c(1,4))

for(i in unique(merged$float)){
  data <- filter(merged, float == i)
  plot(data$chla_argo, -data$pres)
  points(data$t_chla, -data$depth_hplc, col = "green")
  title(i)
  print(i)
}

dev.off()

pdf("plots.pdf", width = 7, height = 7)
d_ply(argo_with_hplc_data, .(float), failwith(NA, function(x){plot(x$chla_adjusted,main=unique(float))}), .print=TRUE)
dev.off()

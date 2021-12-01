library(tidyverse)
library(sf)
library(geosphere)
library(tidybgc)
library(data.table)

#open argo data
argo <- read_csv("Data/argo_db.csv")

#open hplc data
hplc <- read_csv("Data/hplc_argo_campaign")

#remove hplc without position
hplc <- filter(hplc, !is.na(lon) & !is.na(lat))

#correction of a lonlat error around 180 lon
hplc$lat[hplc$lon>100] <- - hplc$lat[hplc$lon>100]
hplc$lon[hplc$lon>100] <- - hplc$lon[hplc$lon>100]

#create a profil id for HPLC data and argo data
hplc <- hplc %>% mutate(profil_id = paste(date, lon, lat, station_name, ctd_number, sep = "_"))
argo <- argo %>% mutate(profil_id = paste(float, lon, lat, sep = "_"))

#create two dataframe with the lon and lat informations
hplc_lonlat <- select(hplc, lon, lat) %>% distinct() %>% na.omit()
argo_lonlat <- select(argo, lon, lat) %>% distinct()

names(hplc_lonlat) <- c("longitude", "latitude")
names(argo_lonlat) <- c("longitude", "latitude")

#empty dataframe to receive information of matched position between argo and hplc data, with the date
match_final <- data.frame("lon_argo" = NA,
                          "lat_argo" = NA,
                          "lon_hplc" = NA,
                          "lat_hplc" = NA)

#for each argo position
for(i in c(1:length(argo_lonlat$longitude))){
  #calculate the distance between the argo profile and all hplc profile
  dist <- distGeo(argo_lonlat[i,c(1,2)], hplc_lonlat)
  #find the minimum distance
  dist_min <- min(dist)
  #find the position of the closest hplc profile
  hplc_position <- hplc_lonlat[which(dist == dist_min),]
  #filter hplc data to extract only the data from the closest profile
  hplc_to_match <- filter(hplc, lon == hplc_position$longitude & lat == hplc_position$latitude)
  #recreate the hplc position df with the date
  hplc_position <- select(hplc_to_match, lon, lat) %>% distinct()
  
  #sometimes we match only one few points on the hplc data
  if(nrow(hplc_to_match) < 3){
    dist_2_min <- order(dist)[1:2]
    #accept profiles that are a bit further than the closest one (1km further) I assume it is due to a drift of the boat during the profile
    hplc_position <- hplc_lonlat[dist_2_min,]
    #extract hplc data to match, using lon, lat and date because of boussole data
    hplc_to_match <- filter(hplc, lon %in% hplc_position$longitude & lat %in% hplc_position$latitude)
    hplc_position <- select(hplc_to_match, lon, lat) %>% distinct()
  }

  #create a dataframe with the argo position and hplc position used to match it
  match <- bind_cols(argo_lonlat[i,], hplc_position)
  names(match) <- c("lon_argo", "lat_argo", "lon_hplc", "lat_hplc")
  #compile it in the final dataframe
  match_final <- bind_rows(match_final, match)
}

#add the position to the closest hplc profile to the argo db
argo_with_hplc_coord_and_date <- argo %>% left_join(match_final, by = c("lon" = "lon_argo", "lat" = "lat_argo"))

#join the hplc data associate to the position added with the line above. And chose the closest depth
argo_with_hplc_data <- argo_with_hplc_coord_and_date %>%
  left_join(hplc, by = c("lon_hplc" = "lon", "lat_hplc" = "lat")) %>% 
  filter(!is.na(chla.x) & !is.na(t_chla)) %>% 
  group_by(lon_hplc, depth.y) %>%
  slice(which.min(abs(pres - depth.y)))

ggplot(argo_with_hplc_data)+
  geom_point(aes(x = t_chla, y = - depth.y))+
  geom_point(aes(x = chla.x, y = - pres, colour = "fluo"))

names(argo_with_hplc_data) <- gsub('.y', '_hplc', colnames(argo_with_hplc_data))
names(argo_with_hplc_data) <- gsub('.x', '_argo', colnames(argo_with_hplc_data))

# Creating new PDF file
pdf("Export_Plots.pdf", width = 16 , height = 10, title = "Argo/HPLC matchup")

# Deciding rows and cols
par(mfrow=c(1,4))

for(i in unique(argo_with_hplc_data$float)){
  data <- filter(argo_with_hplc_data, float == i)
  min_x <- min(c(data$chla_argo, data$t_chla), na.rm = TRUE)
  max_x <- max(c(data$chla_argo, data$t_chla), na.rm = TRUE)
  plot(data$chla_argo, -data$pres, xlim = c(min_x, max_x))
  points(data$t_chla, -data$depth_hplc, col = "green")
  title(i)
  print(i)
}

dev.off()

write_csv(path = "Data/argo_hplc.csv", argo_with_hplc_data)


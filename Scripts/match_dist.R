match_dist <- function(lon_argo, lat_argo, lon_hplc, lat_hplc){
  argo_lonlat <- tibble("lon" = lon_argo, "lat" = lat_argo) %>% distinct()
  hplc_lonlat <- tibble("lon" = lon_hplc, "lat" = lat_hplc) %>% distinct()
  #calculate the distance between the argo profile and all hplc profile
  dist <- distGeo(argo_lonlat, hplc_lonlat)
  
  hplc_lonlat$dist <- dist
  dist <- na.omit(dist)
  
  #find the minimum distance
  dist_min <- min(dist)
  max_dist <- dist_min + 0.1 * dist_min
  min_dist <- dist_min - 0.1 * dist_min
  
  hplc_lonlat <- filter(hplc_lonlat, dist < max_dist & dist > min_dist)
  #find the position of the closest hplc profile
  
  return(hplc_lonlat)
}
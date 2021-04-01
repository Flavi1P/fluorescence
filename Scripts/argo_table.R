library(ncdf4)
library(tidyverse)
library(readxl)
library(lubridate)
library(plyr)
library(ssh)  
library(tidybgc)

#sshfs fpetit@oao2016.obs-vlfr.fr:/data/home oao2016

index_synth <- read_csv("~/oao2016/admt/index_argo/argo_synthetic-profile_index.txt", 
                      skip = 8)

index_bio <- index_synth %>% separate(file, c('dac', 'number', 'prof', 'path'), '/')


ref1 <- read_csv("Data/ref.csv")
refbis <- data.frame("number" = c(6902737, 6902739, 6902735, 6902742, 6902743, 6902880), "lovbio" = c("lovbio103c", "lovbio107c", "lovbio100c", "lovapm002a", "lovapm004a", "lovbio111b"))
refbis$lovbio <- as.character(refbis$lovbio)
ref <- bind_rows(ref1, refbis)
ref <- filter(ref, number != 6901526)

index_bio$number <- as.numeric(index_bio$number)

ref_name <- distinct(left_join(ref, dplyr::select(index_bio, dac, number)))

basepath <- "~/oao2016/admt/GDAC"

ref_name <- ref_name %>%
  mutate(path = paste(basepath, dac, number, 'profiles', paste('SD', number, '_001.nc', sep = ''), sep = '/'))

final_table <- data.frame("depth" = numeric(),
                          "date" = as.Date(x = integer(0), origin = "1970-01-01"),
                          "lon" = numeric(),
                          "lat" = numeric())

for(i in ref_name$path){
  if(file.exists(i)){
    print(paste(i, "exists", sep = " "))
    vars <- index_bio %>% filter(number == str_extract(i, '[0-9]{6,}')) %>% pull(parameters) %>% str_split(" ") %>% unlist() %>% unique()
    table <- extract_sd(i, vars = vars)
  } else{
    print(paste(i, "doesn't exists", sep = " "))
    break()}
  final_table <- bind_rows(final_table, table)
}

final_table <- final_table %>%
  mutate(round_depth = round(pres)) %>%
  group_by(float, round_depth) %>% 
  summarise_all(mean, na.rm = TRUE)

position <- data.frame('lon' = unique(final_table$lon),
                      'lat' = unique(final_table$lat))

worldmap <- map_data('world')

ggplot(position)+
  geom_point(aes(x = lon, y = lat), colour = 'red')+
  geom_polygon(aes(x = long, y = lat, group = group), data = worldmap)+
  coord_quickmap()

write_csv(final_table, 'Data/argo_db.csv')

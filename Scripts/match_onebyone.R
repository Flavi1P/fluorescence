#prototype of a script to match float and hplc manually one by one
library(tidyverse)
library(tidybgc)
library(sf)
library(geosphere)
source("Scripts/read_hplc.R")
source("Scripts/match_dist.R")

#create a list of lovbio to match
#sshfs fpetit@oao2016.obs-vlfr.fr:/data/home oao2016

folders <- list.dirs('.', recursive=TRUE, full.names = TRUE) #list all folders in the raw excel one
HPLC_folders <- folders[grep('PIGMENTS', folders)] #select pigment folder

lovbio_list <- str_extract(HPLC_folders, "/(lov).*|/(tak).*/PIGMENTS")
lovbio_list <- str_extract(lovbio_list, "/.*/")
lovbio_list <- gsub("/", "", lovbio_list)
lovbio_list <- gsub("PIGMENTS", "", lovbio_list)


files <- list.files('./Data/raw_excel', full.names = TRUE)

HPLC_files <- list.files(HPLC_folders, full.names = TRUE) #find path of pigment excel data
HPLC_files_single <- HPLC_files[grep('xls', HPLC_files)]

lov_files <- tibble("lovbio" = lovbio_list, "file" = HPLC_files_single)


#open float database

argo <- read_csv("Data/argo_db.csv")
path_argo <- read_csv("Data/ref_path_argo.csv")

basepath <- "~/oao2016/admt/GDAC"
vars <- c("PRES", "TEMP", "PSAL", "CHLA_ADJUSTED", "BBP700")

ref <- read_csv("Data/ref.csv")
refbis <- data.frame("number" = c(6902737, 6902739, 6902735, 6902742, 6902743, 6902880, 6902734, 6902736, 6902738), "lovbio" = c("lovbio103c", "lovbio107c", "lovbio100c", "lovapm002a", "lovapm004a", "lovbio111b", "lovbio098c", "lovbio101c", "lovbio104c"))

ref <- bind_rows(ref, refbis)
rm(refbis)

path_list <- c()


# lovapm002a --------------------------------------------------------------

#argo file

lovbio_sel <- lov_files$lovbio[1]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_006.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))


#hplc file

hplc_t <- read_hplc(lov_files$file[1])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
    inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "date" = date_argo, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
                gsub("LOD", "0", x)
             }), stringsAsFactors = FALSE)

lovapm002a <- match_finall %>% select(-date, -lovbio)
lovapm002a <- data.frame(lapply(lovapm002a, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovapm002a <- mutate(lovapm002a, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovapm002a)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))



# lovapm004a --------------------------------------------------------------

#argo file

lovbio_sel <- lov_files$lovbio[2]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[2])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovapm004a <- match_finall %>% select(-lovbio)
lovapm004a <- data.frame(lapply(lovapm004a, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovapm004a <- mutate(lovapm004a, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovapm004a)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio001i --------------------------------------------------------------

#argo file

lovbio_sel <- lov_files$lovbio[3]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[3])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- hplc_t %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio001i <- match_finall %>% select(-lovbio)
lovbio001i <- data.frame(lapply(lovbio001i, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio001i <- mutate(lovbio001i, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio001i)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio011b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[5]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[5])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio011b <- match_finall %>% select(-lovbio)
lovbio011b <- data.frame(lapply(lovbio011b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio011b <- mutate(lovbio011b, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio011b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio012b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[6]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[6])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio012b <- match_finall %>% select(-lovbio)
lovbio012b <- data.frame(lapply(lovbio012b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio012b <- mutate(lovbio012b, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio012b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio013b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[7]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[7])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio013b <- match_finall %>% select(-lovbio)
lovbio013b <- data.frame(lapply(lovbio013b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio013b <- mutate(lovbio013b, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio013b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio014b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[8]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[8])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio014b <- match_finall %>% select(-lovbio)
lovbio014b <- data.frame(lapply(lovbio014b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio014b <- mutate(lovbio014b, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio014b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio017b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[9]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[9])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date("2013-04-09")
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio017b <- match_finall %>% select(-lovbio)
lovbio017b <- data.frame(lapply(lovbio017b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio017b <- mutate(lovbio017b, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio017b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio020b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[10]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[10])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio020b <- match_finall %>% select(-lovbio)
lovbio020b <- data.frame(lapply(lovbio020b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio020b <- mutate(lovbio020b, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio020b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio021c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[11]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[11])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio021c <- match_finall %>% select(-lovbio)
lovbio021c <- data.frame(lapply(lovbio021c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio021c <- mutate(lovbio021c, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio021c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio022b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[12]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[12])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio022b <- match_finall %>% select(-lovbio)
lovbio022b <- data.frame(lapply(lovbio022b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio022b <- mutate(lovbio022b, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio022b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio023b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[13]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[13])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio023b <- match_finall %>% select(-lovbio)
lovbio023b <- data.frame(lapply(lovbio023b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE)

lovbio023b <- mutate(lovbio023b, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)

ggplot(lovbio023b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio024c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[14]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[14])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- as.numeric(date_hplc - date_argo)


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio024c <- match_finall %>% select(-lovbio)
lovbio024c <- data.frame(lapply(lovbio024c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio024c <- mutate(lovbio024c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio024c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio025c --------------------------------------------------------------


#argo file

# lovbio_sel <- lov_files$lovbio[15]
# lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)
# 
# path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_007.nc', sep = ''), sep = '/')
# 
# table <- extract_sd(path, vars = vars)
# argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
#   mutate(round_depth = round(pres))
# 
# ggplot(argo_t)+
#   geom_point(aes(x = chla_adjusted, y = - pres))
# #hplc file
# 
# hplc_t <- read_hplc(lov_files$file[13])
# 
# hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)
# 
# dist <- mean(hplc_lonlat$dist)
# 
# hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
# hplc_tomatch$id <- row_number(hplc_tomatch$depth)
# 
# first_match <- argo_t %>%
#   inner_join(hplc_tomatch, by = "round_depth")
# 
# names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
# names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")
# 
# date_argo <- as.Date(first_match$date_argo)
# date_hplc <- as.Date(first_match$date_hplc)
# diff_date <- date_hplc - date_argo
# 
# 
# match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb)
# match_final[match_final == "LOD"] <- "0"
# 
# match_finall <- data.frame(lapply(match_final, function(x) {
#   gsub("LOD", "0", x)
# }), stringsAsFactors = FALSE)
# 
# lovbio023b <- match_finall %>% select(-lovbio)
# lovbio023b <- data.frame(lapply(lovbio023b, function(x) {
#   as.numeric(x)
# }), stringsAsFactors = FALSE)
# 
# lovbio023b <- mutate(lovbio023b, date = date_argo, diff_date = diff_date, lovbio = lovbio_sel, dist = dist)
# 
# ggplot(lovbio023b)+
#   geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
#   geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio026c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[16]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[16])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio026c <- match_finall %>% select(-lovbio)
lovbio026c <- data.frame(lapply(lovbio026c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio026c <- mutate(lovbio026c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio026c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio027b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[17]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[17])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio027b <- match_finall %>% select(-lovbio)
lovbio027b <- data.frame(lapply(lovbio027b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio027b <- mutate(lovbio027b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio027b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio028b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[18]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[18])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio028b <- match_finall %>% select(-lovbio)
lovbio028b <- data.frame(lapply(lovbio028b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio028b <- mutate(lovbio028b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio028b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio029b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[19]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres, 5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[19])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio029b <- match_finall %>% select(-lovbio)
lovbio029b <- data.frame(lapply(lovbio029b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio029b <- mutate(lovbio029b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio029b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio030b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[20]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[20])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio030b <- match_finall %>% select(-lovbio)
lovbio030b <- data.frame(lapply(lovbio030b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio030b <- mutate(lovbio030b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio030b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))



# lovbio031c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[21]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[21])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio031c <- match_finall %>% select(-lovbio)
lovbio031c <- data.frame(lapply(lovbio031c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio031c <- mutate(lovbio031c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio031c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio032b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[22]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[22])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio032b <- match_finall %>% select(-lovbio)
lovbio032b <- data.frame(lapply(lovbio032b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio032b <- mutate(lovbio032b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio032b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio035b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[23]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[23])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio035b <- match_finall %>% select(-lovbio)
lovbio035b <- data.frame(lapply(lovbio035b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio035b <- mutate(lovbio035b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio035b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio037c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[24]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[24])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio037c <- match_finall %>% select(-lovbio)
lovbio037c <- data.frame(lapply(lovbio037c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio037c <- mutate(lovbio037c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio037c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio038b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[25]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[25])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio038b <- match_finall %>% select(-lovbio)
lovbio038b <- data.frame(lapply(lovbio038b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio038b <- mutate(lovbio038b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio038b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio040b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[26]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres, 5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[26])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio040b <- match_finall %>% select(-lovbio)
lovbio040b <- data.frame(lapply(lovbio040b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio040b <- mutate(lovbio040b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio040b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio042d --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[27]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_003.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[27])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio042d <- match_finall %>% select(-lovbio)
lovbio042d <- data.frame(lapply(lovbio042d, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio042d <- mutate(lovbio042d, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio042d)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio043b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[28]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_003.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[28])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio043b <- match_finall %>% select(-lovbio)
lovbio043b <- data.frame(lapply(lovbio043b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio043b <- mutate(lovbio043b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio043b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio044b --------------------------------------------------------------


#argo file
# 
# lovbio_sel <- lov_files$lovbio[29]
# lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)
# 
# path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_013.nc', sep = ''), sep = '/')
# 
# table <- extract_sd(path, vars = vars)
# argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
#   mutate(round_depth = round(pres))
# 
# ggplot(argo_t)+
#   geom_point(aes(x = chla_adjusted, y = - pres))
# #hplc file
# 
# hplc_t <- read_hplc(lov_files$file[28])
# 
# hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)
# 
# dist <- mean(hplc_lonlat$dist)
# 
# hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
# hplc_tomatch$id <- row_number(hplc_tomatch$depth)
# 
# first_match <- argo_t %>%
#   inner_join(hplc_tomatch, by = "round_depth")
# 
# names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
# names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")
# 
# date_argo <- as.Date(first_match$date_argo)
# date_hplc <- as.Date(first_match$date_hplc)
# diff_date <- date_hplc - date_argo
# 
# 
# match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
# match_final[match_final == "LOD"] <- "0"
# 
# match_finall <- data.frame(lapply(match_final, function(x) {
#   gsub("LOD", "0", x)
# }), stringsAsFactors = FALSE)
# 
# lovbio043b <- match_finall %>% select(-lovbio)
# lovbio043b <- data.frame(lapply(lovbio043b, function(x) {
#   as.numeric(x)
# }), stringsAsFactors = FALSE) %>% 
#   group_by(round_depth) %>% 
#   summarise_all(mean) %>% 
#   ungroup() %>% 
#   select(-round_depth)
# 
# lovbio043b <- mutate(lovbio043b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)
# 
# ggplot(lovbio043b)+
#   geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
#   geom_point(aes(x = t_chla, y = - depth, colour = "chla")) #no chla adjusted
# lovbio045b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[30]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres, 5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[30])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio045b <- match_finall %>% select(-lovbio)
lovbio045b <- data.frame(lapply(lovbio045b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio045b <- mutate(lovbio045b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio045b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio050b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[31]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[31])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio050b <- match_finall %>% select(-lovbio)
lovbio050b <- data.frame(lapply(lovbio050b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio050b <- mutate(lovbio050b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio050b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio057b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[32]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[32])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio057b <- match_finall %>% select(-lovbio)
lovbio057b <- data.frame(lapply(lovbio057b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio057b <- mutate(lovbio057b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio057b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio059c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[33]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[33])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date("2014-06-15")
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio059c <- match_finall %>% select(-lovbio)
lovbio059c <- data.frame(lapply(lovbio059c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio059c <- mutate(lovbio059c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio059c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio061c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[34]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[34])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date("2014-06-09")
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio061c <- match_finall %>% select(-lovbio)
lovbio061c <- data.frame(lapply(lovbio061c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio061c <- mutate(lovbio061c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio061c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio063c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[35]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres, 5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[35])

hplc_t <- filter(hplc_t, date == "2014-07-07")

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio063c <- match_finall %>% select(-lovbio)
lovbio063c <- data.frame(lapply(lovbio063c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio063c <- mutate(lovbio063c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio063c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio064b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[36]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[36])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio064b <- match_finall %>% select(-lovbio)
lovbio064b <- data.frame(lapply(lovbio064b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio064b <- mutate(lovbio064b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio064b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio067c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[38]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[38])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio067c <- match_finall %>% select(-lovbio)
lovbio067c <- data.frame(lapply(lovbio067c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio067c <- mutate(lovbio067c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio067c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio068d --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[39]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[39])

hplc_t <- filter(hplc_t, date == "2014-07-07")

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio068d <- match_finall %>% select(-lovbio)
lovbio068d <- data.frame(lapply(lovbio068d, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio068d <- mutate(lovbio068d, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio068d)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio075b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[40]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[40])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio075b <- match_finall %>% select(-lovbio)
lovbio075b <- data.frame(lapply(lovbio075b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio075b <- mutate(lovbio075b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio075b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio077b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[41]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[41])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio077b <- match_finall %>% select(-lovbio)
lovbio077b <- data.frame(lapply(lovbio077b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio077b <- mutate(lovbio077b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio077b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))



# lovbio079b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[42]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[42])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio079b <- match_finall %>% select(-lovbio)
lovbio079b <- data.frame(lapply(lovbio079b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio079b <- mutate(lovbio079b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio079b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))




# lovbio082b --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[43]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,2))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[43])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio082b <- match_finall %>% select(-lovbio)
lovbio082b <- data.frame(lapply(lovbio082b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio082b <- mutate(lovbio082b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio082b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio083d --------------------------------------------------------------


# #argo file
# 
# lovbio_sel <- lov_files$lovbio[44]
# lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)
# 
# path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
# 
# table <- extract_sd(path, vars = vars)
# argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
#   mutate(round_depth = round(pres))
# 
# ggplot(argo_t)+
#   geom_point(aes(x = chla_adjusted, y = - pres))+
#   ylim(-50,0)
# #hplc file
# 
# hplc_t <- read_hplc(lov_files$file[44])
# 
# hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)
# 
# dist <- mean(hplc_lonlat$dist)
# 
# hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
# hplc_tomatch$id <- row_number(hplc_tomatch$depth)
# 
# first_match <- argo_t %>%
#   inner_join(hplc_tomatch, by = "round_depth")
# 
# names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
# names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")
# 
# date_argo <- as.Date(first_match$date_argo)
# date_hplc <- as.Date(first_match$date_hplc)
# diff_date <- date_hplc - date_argo
# 
# 
# match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
# match_final[match_final == "LOD"] <- "0"
# 
# match_finall <- data.frame(lapply(match_final, function(x) {
#   gsub("LOD", "0", x)
# }), stringsAsFactors = FALSE)
# 
# lovbio077b <- match_finall %>% select(-lovbio)
# lovbio077b <- data.frame(lapply(lovbio077b, function(x) {
#   as.numeric(x)
# }), stringsAsFactors = FALSE) %>% 
#   group_by(round_depth) %>% 
#   summarise_all(mean) %>% 
#   ungroup() %>% 
#   select(-round_depth)
# 
# lovbio077b <- mutate(lovbio077b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)
# 
# ggplot(lovbio077b)+
#   geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
#   geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


 #no match with xcel bioargomed
# lovbio084d --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[45]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[45])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio084d <- match_finall %>% select(-lovbio)
lovbio084d <- data.frame(lapply(lovbio084d, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio084d <- mutate(lovbio084d, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio084d)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))
 #DCM shift on the two profiles
# lovbio085d --------------------------------------------------------------


#argo file
# 
# lovbio_sel <- lov_files$lovbio[46]
# lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)
# 
# path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
# 
# table <- extract_sd(path, vars = vars)
# argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
#   mutate(round_depth = plyr::round_any(pres,2))
# 
# ggplot(argo_t)+
#   geom_point(aes(x = chla_adjusted, y = - pres))
# #hplc file
# 
# hplc_t <- read_hplc(lov_files$file[46])
# 
# hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)
# 
# dist <- mean(hplc_lonlat$dist)
# 
# hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
# hplc_tomatch$id <- row_number(hplc_tomatch$depth)
# 
# first_match <- argo_t %>%
#   inner_join(hplc_tomatch, by = "round_depth")
# 
# names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
# names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")
# 
# date_argo <- as.Date(first_match$date_argo)
# date_hplc <- as.Date(first_match$date_hplc)
# diff_date <- date_hplc - date_argo
# 
# 
# match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
# match_final[match_final == "LOD"] <- "0"
# 
# match_finall <- data.frame(lapply(match_final, function(x) {
#   gsub("LOD", "0", x)
# }), stringsAsFactors = FALSE)
# 
# lovbio082b <- match_finall %>% select(-lovbio)
# lovbio082b <- data.frame(lapply(lovbio082b, function(x) {
#   as.numeric(x)
# }), stringsAsFactors = FALSE) %>% 
#   group_by(round_depth) %>% 
#   summarise_all(mean) %>% 
#   ungroup() %>% 
#   select(-round_depth)
# 
# lovbio082b <- mutate(lovbio082b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)
# 
# ggplot(lovbio082b)+
#   geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
#   geom_point(aes(x = t_chla, y = - depth, colour = "chla")) #no match with xcel bioargomed
# lovbio086d --------------------------------------------------------------


# #argo file
# 
# lovbio_sel <- lov_files$lovbio[47]
# lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)
# 
# path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_040.nc', sep = ''), sep = '/')
# 
# table <- extract_sd(path, vars = vars)
# argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
#   mutate(round_depth = plyr::round_any(pres,2))
# 
# ggplot(argo_t)+
#   geom_point(aes(x = chla_adjusted, y = - pres))
# #hplc file
# 
# hplc_t <- read_hplc(lov_files$file[43])
# 
# hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)
# 
# dist <- mean(hplc_lonlat$dist)
# 
# hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
# hplc_tomatch$id <- row_number(hplc_tomatch$depth)
# 
# first_match <- argo_t %>%
#   inner_join(hplc_tomatch, by = "round_depth")
# 
# names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
# names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")
# 
# date_argo <- as.Date(first_match$date_argo)
# date_hplc <- as.Date(first_match$date_hplc)
# diff_date <- date_hplc - date_argo
# 
# 
# match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
# match_final[match_final == "LOD"] <- "0"
# 
# match_finall <- data.frame(lapply(match_final, function(x) {
#   gsub("LOD", "0", x)
# }), stringsAsFactors = FALSE)
# 
# lovbio082b <- match_finall %>% select(-lovbio)
# lovbio082b <- data.frame(lapply(lovbio082b, function(x) {
#   as.numeric(x)
# }), stringsAsFactors = FALSE) %>% 
#   group_by(round_depth) %>% 
#   summarise_all(mean) %>% 
#   ungroup() %>% 
#   select(-round_depth)
# 
# lovbio082b <- mutate(lovbio082b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)
# 
# ggplot(lovbio082b)+
#   geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
#   geom_point(aes(x = t_chla, y = - depth, colour = "chla")) #no chla adjusted
# lovbio088d --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[48]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_004.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = round(pres))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[48])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio088d <- match_finall %>% select(-lovbio)
lovbio088d <- data.frame(lapply(lovbio088d, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio088d <- mutate(lovbio088d, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio088d)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio089d --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[49]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,2))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[49])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio089d <- match_finall %>% select(-lovbio)
lovbio089d <- data.frame(lapply(lovbio089d, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio089d <- mutate(lovbio089d, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio089d)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio090d --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[50]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,2))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[50])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= round(depth), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio082b <- match_finall %>% select(-lovbio)
lovbio082b <- data.frame(lapply(lovbio082b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio082b <- mutate(lovbio082b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio082b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla")) # bioargomed pourri
# lovbio091d --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[51]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,2))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[51])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio091b <- match_finall %>% select(-lovbio)
lovbio091b <- data.frame(lapply(lovbio091b, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio091b <- mutate(lovbio091b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio091b)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio093d --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[52]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,2))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[52])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio093d <- match_finall %>% select(-lovbio)
lovbio093d <- data.frame(lapply(lovbio093d, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio093d <- mutate(lovbio093d, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio093d)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio098c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[53]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,2))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[53])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio098c <- match_finall %>% select(-lovbio)
lovbio098c <- data.frame(lapply(lovbio098c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio098c <- mutate(lovbio098c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio098c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))


# lovbio100c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[54]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[54])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio100c <- match_finall %>% select(-lovbio)
lovbio100c <- data.frame(lapply(lovbio100c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio100c <- mutate(lovbio100c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio100c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio101c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[55]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[55])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio101c <- match_finall %>% select(-lovbio)
lovbio101c <- data.frame(lapply(lovbio101c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio101c <- mutate(lovbio101c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio101c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio103c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[56]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[56])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio103c <- match_finall %>% select(-lovbio)
lovbio103c <- data.frame(lapply(lovbio103c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio103c <- mutate(lovbio103c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio103c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio104c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[58]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[58])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio104c <- match_finall %>% select(-lovbio)
lovbio104c <- data.frame(lapply(lovbio104c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio104c <- mutate(lovbio104c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio104c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio107c --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[59]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_002.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_hplc(lov_files$file[59])

hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)

dist <- mean(hplc_lonlat$dist)

hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
hplc_tomatch$id <- row_number(hplc_tomatch$depth)

first_match <- argo_t %>%
  inner_join(hplc_tomatch, by = "round_depth")

names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")

date_argo <- as.Date(first_match$date_argo)
date_hplc <- as.Date(first_match$date_hplc)
diff_date <- date_hplc - date_argo


match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
match_final[match_final == "LOD"] <- "0"

match_finall <- data.frame(lapply(match_final, function(x) {
  gsub("LOD", "0", x)
}), stringsAsFactors = FALSE)

lovbio107c <- match_finall %>% select(-lovbio)
lovbio107c <- data.frame(lapply(lovbio107c, function(x) {
  as.numeric(x)
}), stringsAsFactors = FALSE) %>% 
  group_by(round_depth) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  select(-round_depth)

lovbio107c <- mutate(lovbio107c, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)

ggplot(lovbio107c)+
  geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "chla"))

# lovbio111b --------------------------------------------------------------

# 
# #argo file
# 
# lovbio_sel <- lov_files$lovbio[60]
# lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)
# 
# path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
# 
# table <- extract_sd(path, vars = vars)
# argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
#   mutate(round_depth = plyr::round_any(pres,5))
# 
# ggplot(argo_t)+
#   geom_point(aes(x = chla_adjusted, y = - pres))
# #hplc file
# 
# hplc_t <- read_hplc(lov_files$file[60])
# 
# hplc_lonlat <- match_dist(argo_t$lon, argo_t$lat, hplc_t$lon, hplc_t$lat)
# 
# dist <- mean(hplc_lonlat$dist)
# 
# hplc_tomatch <- filter(hplc_t, lon %in% hplc_lonlat$lon, lat %in% hplc_lonlat$lat) %>% mutate(round_depth= plyr::round_any(depth, 5), lovbio = lovbio_sel)
# hplc_tomatch$id <- row_number(hplc_tomatch$depth)
# 
# first_match <- argo_t %>%
#   inner_join(hplc_tomatch, by = "round_depth")
# 
# names(first_match) <- str_replace(colnames(first_match), "\\.x", "_argo")
# names(first_match) <- str_replace(colnames(first_match), "\\.y", "_hplc")
# 
# date_argo <- as.Date(first_match$date_argo)
# date_hplc <- as.Date(first_match$date_hplc)
# diff_date <- date_hplc - date_argo
# 
# 
# match_final <- first_match %>% select(lovbio, "float_num" = float, "lon" = lon_argo, "lat" = lat_argo, "depth" = pres, "fluo" = chla_adjusted, chla, t_chla, peri, but, hex, zea, fuco, allo, chlb, round_depth)
# match_final[match_final == "LOD"] <- "0"
# 
# match_finall <- data.frame(lapply(match_final, function(x) {
#   gsub("LOD", "0", x)
# }), stringsAsFactors = FALSE)
# 
# lovbio111b <- match_finall %>% select(-lovbio)
# lovbio111b <- data.frame(lapply(lovbio111b, function(x) {
#   as.numeric(x)
# }), stringsAsFactors = FALSE) %>% 
#   group_by(round_depth) %>% 
#   summarise_all(mean) %>% 
#   ungroup() %>% 
#   select(-round_depth)
# 
# lovbio111b <- mutate(lovbio111b, date = unique(date_argo), diff_date = unique(diff_date), lovbio = lovbio_sel, dist = dist)
# 
# ggplot(lovbio111b)+
#   geom_point(aes(x = fluo, y = - depth, colour = "fluo"))+
#   geom_point(aes(x = t_chla, y = - depth, colour = "chla"))
 #unmatch between hplc sheet and argo file (soclim)
# takuvik --------------------------------------------------------------


#argo file

lovbio_sel <- lov_files$lovbio[61]
lovbio_path <-  filter(path_argo, lovbio == lovbio_sel)

path = paste(basepath, lovbio_path$dac, lovbio_path$number, 'profiles', paste('SD', lovbio_path$number, '_001.nc', sep = ''), sep = '/')
path_list <- c(path_list, path)

table <- extract_sd(path, vars = vars)
argo_t <- filter(table, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8)) %>% 
  mutate(round_depth = plyr::round_any(pres,5))

ggplot(argo_t)+
  geom_point(aes(x = chla_adjusted, y = - pres))
#hplc file

hplc_t <- read_excel(lov_files$file[61]) %>% janitor::clean_names() %>% 
  mutate(date = sampling_date_utc,
         lovbio = x11,
         lon = longitude,
         lat = latitude,
         t_chla = total_chlorophyll_a_hplc_mg_m_3,
         fluo = chla_adjusted) %>% 
  select(date, lovbio, lon, lat, depth, t_chla, fluo)

takuvik <- hplc_t %>% mutate(peri = NA,
                             chla = NA,
                             but = NA,
                             hex = NA,
                             zea = NA,
                             fuco = NA,
                             allo = NA,
                             chlb = NA,
                             diff_date = 0,
                             dist = NA, 
                             depth = - depth)
takuvik <- inner_join(takuvik, ref) %>% dplyr::rename("float_num" = number)


# bind dataframe ----------------------------------------------------------


argo_matchup <- lovapm002a %>% mutate(diff_date = as.numeric(diff_date))
for(i in list(lovapm002a, lovapm004a, lovbio001i, lovbio011b, lovbio012b, lovbio013b, lovbio014b, lovbio017b, lovbio020b,
              lovbio021c, lovbio022b, lovbio023b, lovbio024c, lovbio026c, lovbio027b, lovbio028b, lovbio029b, lovbio030b,
              lovbio031c, lovbio032b, lovbio035b, lovbio037c, lovbio038b, lovbio040b, lovbio042d, lovbio043b, lovbio045b,
              lovbio050b, lovbio057b, lovbio059c, lovbio061c, lovbio063c, lovbio064b, lovbio064b, lovbio067c, lovbio068d, lovbio075b,
              lovbio077b, lovbio079b, lovbio082b, lovbio084d, lovbio088d, lovbio089d,
              lovbio091b, lovbio093d, lovbio098c, lovbio100c, lovbio101c, lovbio103c, lovbio104c, lovbio107c,
              takuvik)){
  t <- i %>% mutate(diff_date = as.numeric(diff_date))
  argo_matchup <- bind_rows(argo_matchup, t)
}

ggplot(argo_matchup)+
  geom_point(aes(x = fluo * 2, y = -depth, colour = "fluo"))+
  geom_point(aes(x = t_chla, y = - depth, colour = "hplc"))+
  facet_wrap(.~lovbio, scales = "free")

worldmap <- map_data("world")

ggplot(argo_matchup)+
  geom_point(aes(x = lon, y = lat, colour = lovbio))+
  geom_polygon(aes(x = long, y = lat, group = group), data = worldmap)+
  coord_quickmap()

#write_csv(argo_matchup, "Data/argo_matchup.csv")

argo_db <- tibble()

for(i in path_list){
  argo_t <- extract_sd(i, vars = vars)
  argo_t <- filter(argo_t, pres_qc %in% c(1,2,5,8), temp_qc %in% c(1,2,5,8), psal_qc %in% c(1,2,5,8), chla_adjusted_qc %in% c(1,2,5,8))
  argo_db <- bind_rows(argo_db, argo_t)
}

argo_db$float <- as.numeric(argo_db$float)
argo_dbb <- left_join(argo_db, argo_matchup, by = c("float" = "float_num"))

ggplot(argo_db)+
  geom_point(aes(x = chla_adjusted, y = -pres, colour = "matchup"))+
  facet_wrap(.~float, scales = "free")+
  ylim(-200, 0)

source("Scripts/zeu_moma.R")

ze_vec <- c()
for(i in unique(argo_db$float)){
  t <- filter(argo_db, float == i)
  ze <- Zeu_moma(t$chla_adjusted, t$pres)
  ze_vec <- c(ze_vec, ze)
}

ze_table <- data.frame("float" = unique(argo_db$float), "ze" = ze_vec)

argo_dbb <- left_join(argo_db, ze_table)

ggplot(argo_dbb)+
  geom_point(aes(x = chla_adjusted, y = -pres, colour = "matchup"))+
  geom_hline(aes(yintercept = -ze))+
  facet_wrap(.~float, scales = "free")+
  ylim(-200, 0)

argo_matchup_with_ze <- left_join(argo_matchup, ze_table, by = c("float_num" = "float"))

table(is.na(argo_matchup_with_ze$ze))

#write_csv(argo_matchup_with_ze, "Data/argo_matchup_ze.csv")

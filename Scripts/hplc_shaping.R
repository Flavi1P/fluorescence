library(readxl)
library(tidyverse, lubridate)

#extract path to HPLC excel sheets, only one sheet per campagne ####
folders <- list.dirs('.', recursive=TRUE, full.names = TRUE) #list all folders in the raw excel one
HPLC_folders <- folders[grep('PIGMENTS', folders)] #select pigment folder
files <- list.files('./Data/raw_excel', full.names = TRUE)
HPLC_files_single <- files[grep('xls', files)]
HPLC_files <- list.files(HPLC_folders, full.names = TRUE) #find path of pigment excel data
unique_hplc <- unique(gsub('\\./.+/', '', HPLC_files))#extract unique file names

my_path <- c()
check <- c()
for(i in HPLC_files){
  test <- gsub('\\./.+/', '', i)
  if (test %in% check){}
  else{
    my_path <- append(my_path, i)
    check <- append(check, test)
  }
} #remove duplicate files from my paths

my_path <- my_path[grep('\\.xls', my_path)]
my_path <- my_path[-grep('autre', my_path)]
my_path <- c(my_path, HPLC_files_single)
rm(files, HPLC_files_single, folder, HPLC_folders, HPLC_files, unique_hplc, check, i, test, folders)

#create a df with all pigments ####
cruises_data <- tibble('path' = my_path,
                       'cruise' = c('soclim', 'boussole', 'labrador', 'icb', 'dewex', 'irs', 'dewex', 'oiso', 'lovbio059c', 'lovbio061c', 'moose', 'outpace', 'bioargomed', 'liste_flotteurs', 'greenedge', 'mobydick', 'soclim'))

data_large <- tibble('cruise_name' = character(),
                   'date' = character(),
                   'lat' = numeric(),
                   'lon' = numeric(),
                   'station_name' = character(),
                   'ctd_number' = character(),
                   'bottle_number' = numeric(),
                   'depth' = numeric(),
                   'chl_c1_c2' = character(),
                    'chl_c3' = character(),
               'peri' = character(),
               'but' = character(),
               'fuco' = character(),
               'neox' = character(),
               'prasi' = character(),
               'viola' = character(),
               'hex' = character(),
               'diad' = character(),
               'allo' = character(),
               'diat' = character(),
               'zea' = character(),
               'lutein' = character(),
               'dv_chlb' = character(),
               'chlb' = character(),
               'tchlb' = character(),
               'dv_chla' = character(),
               'chla' = character(),
               'tchla' = character(),
               'phaeo' = character())
for (i in 1:length(cruises_data$path)){
  path <- as.character(cruises_data[i,1])
  
  data <- read_excel(path) %>% janitor::clean_names()

  cruise_name <- cruises_data[i,2]
  
  date <- select(data, matches('(^.+date)|(^date$)'))
  if(length(date) != 1){
    if(length(date) == 2){
      date <- select(date, matches('utc'))
    }
    if(length(date) == 0){
      year <- select(data, matches('year'))
      month <- select(data, matches('month'))
      day <- select(data, matches('day'))
      date <- paste(day[[1]], month[[1]], year[[1]], sep = '/')
    }
  }
  
  lat <- select(data, matches('lat'))
  if(length(lat)!=1){
    lat <- lat[,1]
    if(length(lat) != 1){
      print(i)
      break()
    }
  }
  
  lon <- select(data, matches('lon'))
  if(length(lon) != 1){
    lon <- lon[,1]
    if(length(lon) != 1){
      print(i)
      break()
    }
  }
  
  ctd <- select(data, matches('ctd'))
  if(length(ctd) != 1){
    print(i)
    break()
  }
  
  station <- select(data, matches('station'))
  if(length(station) != 1){
    station <- date
    print(paste('station = date for ', path))
  }
  
  btl <- select(data, matches('bottle'))
  if(length(btl) != 1){
    print(paste('bottle is NA for ', path))
    btl <- NA
  }  
  
  depth <- select(data, matches('depth'))
  if(length(depth) != 1){
    depth <- select(data, matches('^depth$'))
    if(length(depth) != 1){
      print(path)
      break()
    }

  }
  chl_c1_c2 <- NA
  chl_c3 <- NA
  peri <- NA
  but <- NA
  fuco <- NA
  neox <- NA
  prasi <- NA
  viola <- NA
  hex <- NA
  diad <- NA
  allo <- NA
  diat <- NA
  zea <- NA
  lutein <- NA
  dv_chlb <- NA
  chlb <- NA
  tchlb <- NA
  dv_chla <- NA
  chla <- NA
  tchla <- NA
  if(i != 14){
    
    chl_c1_c2 <- select(data, matches('(chl.*c1.*c2)$|(chl.*c2.*c1)$|(chl.*c2.*c1).*[^a]$|(chl.*c1.*c2).*[^a]$'))
    if(length(chl_c1_c2) != 1){
      chl_c1 <- select(data, matches('chl.+c1$'))
      chl_c2 <- select(data, matches('chl.+c2$'))
      chl_c1_c2 <- chl_c1 + chl_c2
      if(length(chl_c1) !=1 | length(chl_c2) !=1){
        print(path)
        break() 
      }
    }
    
    chl_c3 <- select(data, matches('^chl.*c3$'))
    if(length(chl_c3) !=1){
      print(paste('wrong regex for chl c3 on ', path))
      break()
    }
    
    peri <- select(data, matches('peridinin'))
    if(length(peri) !=1){
      print(paste('wrong regex for peri on ', path))
      break()
    }
    
    but <- select(data, matches('(?=.*but)(?!.*qa)', perl = TRUE))
    if(length(but) !=1){
      print(paste('wrong regex for but on ', path))
      break()
    }
    
    fuco <- select(data, matches('^(?=fuco)(?!.*qa)', perl = TRUE))
    if(length(but) !=1){
      print(paste('wrong regex for fuco on ', path))
      break()
    } 
    
    neox <- select(data, matches('^(?=neox)(?!.*qa)', perl = TRUE))
    if(length(neox) !=1){
      print(paste('wrong regex for neox on ', path))
      break()
    }
    
    prasi <- select(data, matches('^(?=prasi)(?!.*qa)', perl = TRUE))
    if(length(prasi) !=1){
      print(paste('wrong regex for prasi on ', path))
      break()
    }
    
    viola <- select(data, matches('^(?=viola)(?!.*qa)', perl = TRUE))
    if(length(viola) !=1){
      print(paste('wrong regex for viola on ', path))
      break()
    }
    
    hex <- select(data, matches('(?=hex)(?!.*qa)', perl = TRUE))
    if(length(hex) !=1){
      print(paste('wrong regex for hex on ', path))
      break()
    }
    
    diad <- select(data, matches('(?=diad)(?!.*qa)', perl = TRUE))
    if(length(diad) !=1){
      print(paste('wrong regex for diad on ', path))
      break()
    }
    
    allo <- select(data, matches('(?=allo)(?!.*qa)', perl = TRUE))
    if(length(allo) !=1){
      print(paste('wrong regex for allo on ', path))
      break()
    }
    
    diat <- select(data, matches('(?=diat)(?!.*qa)', perl = TRUE))
    if(length(diat) !=1){
      print(paste('wrong regex for diat on ', path))
      break()
    }
    
    zea <- select(data, matches('(?=zea)(?!.*qa)', perl = TRUE))
    if(length(allo) !=1){
      print(paste('wrong regex for zea on ', path))
      break()
    }
    
    lutein <- select(data, matches('^(?=lutein)(?!.*qa)', perl = TRUE))
    if(length(lutein) !=1){
      print(paste('wrong regex for lutein on ', path))
      break()
    }
    
    dv_chlb <-  select(data, matches('(?=div.*b$)(?!.*qa)', perl = TRUE))
    if(length(dv_chlb) !=1){
      if(cruise_name %in% c('dewex', 'lovbio059c', 'lovbio061c', 'lovbio063c', 'moose', 'greenedge')){
        dv_chlb <- NA
      }
      else{
        print(paste('wrong regex for dv_chlb on ', path))
        break()
      }
    }
    
    chlb <- select(data, matches('^(?=chl.*b$)(?!.*qa)', perl = TRUE))
    if(length(chlb) !=1){
      if(cruise_name %in% c('dewex', 'lovbio059c', 'lovbio061c', 'lovbio063c', 'moose')){
        chlb <- NA
      }
      else{
      print(paste('wrong regex for chlb on ', path))
      break()
      }
    }
    
    tchlb <- select(data, matches('^(?=t.*chl.*b$)(?!.*qa)', perl = TRUE))
    if(cruise_name == 'greenedge'){
      tchlb <- chlb
    }
    if(length(tchlb) !=1){
      print(paste('wrong regex for tchlb on ', path))
      break()
    }
    
    dv_chla <- select(data, matches('(?=div.*a$)(?!.*qa)', perl = TRUE))
    if(length(dv_chla) !=1){
      print(paste('wrong regex for dv_chla on ', path))
      break()
    }
    
    chla <- select(data, matches('^(?!chlorophyllid)(?=chl.*a$)(?!.*qa)', perl = TRUE))
    if(length(chla) !=1){
      print(paste('wrong regex for chla on ', path))
      break()
    }
    
    tchla <- select(data, matches('^(?=t.*chl.*a$)(?!.*qa)', perl = TRUE))
    if(length(tchla) !=1){
      print(paste('wrong regex for tchla on ', path))
      break()
    }
    
    phaeophytin <- select(data, matches('(?=phaeophytin)(?!.*qa)', perl = TRUE))
    if(length(phaeophytin) !=1){
      print(paste('wrong regex for phaeophytin on ', path))
      break()
    }
    
  }
  
  
  data_interm <- tibble('cruise_name' = cruise_name[[1]],
                        'date' = as.character(date[[1]]),
                        'lon' = as.numeric(lon[[1]]),
                        'lat' = as.numeric(lat[[1]]),
                        'ctd_number' = as.character(ctd[[1]]),
                        'station_name' = as.character(station[[1]]),
                        'bottle_number' = as.numeric(btl[[1]]),
                        'chl_c1_c2' = as.character(chl_c1_c2[[1]]),
                        'chl_c3' = as.character(chl_c3[[1]]),
                        'peri' = as.character(peri[[1]]),
                        'but' = as.character(but[[1]]),
                        'fuco' = as.character(fuco[[1]]),
                        'neox' = as.character(neox[[1]]),
                        'prasi' = as.character(prasi[[1]]),
                        'viola' = as.character(viola[[1]]),
                        'hex' = as.character(hex[[1]]),
                        'diad' = as.character(diad[[1]]),
                        'allo' = as.character(allo[[1]]),
                        'diat' = as.character(diat[[1]]),
                        'zea' = as.character(zea[[1]]),
                        'lutein' = as.character(lutein[[1]]),
                        'dv_chlb' = as.character(dv_chlb[[1]]),
                        'chlb' = as.character(chlb[[1]]),
                        't_chlb' = as.character(tchlb[[1]]),
                        'dv_chla' = as.character(dv_chla[[1]]),
                        'chla' = as.character(chla[[1]]),
                        't_chla' = as.character(tchla[[1]]),
                        'depth' = abs(as.numeric(depth[[1]])))
  data_large <- bind_rows(data_large, data_interm)
  
}

data_long <- data_large %>% pivot_longer(chl_c1_c2:t_chla, names_to = 'pigment', values_to = 'concentration')

data_long$concentration <- tolower(data_long$concentration)
data_long$concentration[data_long$concentration == 'lod'] <- 0



table(is.na(data_long$concentration))
data_long$concentration <- as.numeric(data_long$concentration)

data_long$concentration[data_long$concentration < 0] <- NA



date_num <- data_large$date[grep('^[0-9]{5}$', data_large$date)]


dates <- unique(data_large$date)

date_ok <-  data_large$date[grep('^([0-9]{4}\\-(0?[1-9]|1[0-2])\\-(0?[1-9]|[1-2][0-9]|3[0-1]))$', data_large$date)]

date_ok <- as_date(date_ok)
date_ok <- unique(date_ok)

date_num_date <- as.Date(as.numeric(date_num), origin = '1899-10-31')
data_large$date[grep('^[0-9]{5}$', data_large$date)] <- as.character(date_num_date)



other_format <- unique(data_large$date[-grep('^([0-9]{4}\\-(0?[1-9]|1[0-2])\\-(0?[1-9]|[1-2][0-9]|3[0-1]))$', data_large$date)])
convert_df <- data.frame('good_date' = c("2001-07-21", '2013-02-13', '2014-06-15', '2014-06-09'),
                         'bad_date' = c("24/7/2001", "13/0/2013", "15/06/2014", "09/06/2014"))
other_format <- (data_large$date[-grep('^([0-9]{4}\\-(0?[1-9]|1[0-2])\\-(0?[1-9]|[1-2][0-9]|3[0-1]))$', data_large$date)])

df_to_convert <- data.frame('bad_date' = other_format)
df_date <- left_join(df_to_convert, convert_df)

data_large$date[-grep('^([0-9]{4}\\-(0?[1-9]|1[0-2])\\-(0?[1-9]|[1-2][0-9]|3[0-1]))$', data_large$date)] <- as.character(df_date$good_date)

data_large$date <- as.Date(data_large$date)

ggplot(data_large)+
  geom_density(aes(x = date))

world_map <- map_data('world')

ggplot()+
  geom_polygon(aes(x = long, y = lat, group = group), data = world_map)+
  geom_point(aes(x = lon, y = lat, colour = cruise_name), data = data_large)+
  coord_quickmap()


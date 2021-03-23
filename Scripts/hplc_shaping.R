library(readxl)
library(tidyverse)

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
my_path <- c(my_path, HPLC_files_single)
rm(files, HPLC_files_single, folder, HPLC_folders, HPLC_files, unique_hplc, check, i, test, folders)

#store all possibles colnames ####
cruises_data <- tibble('path' = my_path,
                       'cruise' = c('soclim', 'lovbio001i', 'boussole', 'labrador', 'icb', 'dewex', 'irs', 'dewex', 'oiso', 'lovbio059c', 'lovbio061c', 'moose', 'lovbio064b', 'outpace', 'bioargomed', 'liste_flotteurs', 'greenedge', 'mobydick', 'soclim'))

list_colname <- c()
test <- tibble('cruise_name' = character(),
                   'date' = character(),
                   'lat' = numeric(),
                   'lon' = numeric(),
                   'station_name' = character(),
                   'ctd_number' = numeric(),
                   'bottle_number' = numeric(),
                   'depth' = numeric(),
                   'pigment' = character(),
                   'concentration' = character(),
                   'flag' = character()
                   )
for (i in length(cruises_data$path)){
  path <- as.character(cruises_data[i,1])
  
  data <- read_excel(path) %>% janitor::clean_names()
  # column <- colnames(data)
  # list_colname <- append(list_colname, column)
  # list_colname <- unique(list_colname)
  cruise_name <- cruises_data[i,2]
  date <- select(data, matches('(^.+date)|date'))
  ship <- select(data, matches('r/v|ship'))
  
  
  data_to_bind <- tibble('cruise_name' = cruise_name)
  
  test <- bind_rows(test, data)
}

check_colname <- function(regexp, data){
  result <- data[grep(regexp, data)]
  return(result)
}



check_colname('^t.*[^q]a$', list_colname)
check_colname('^(chl)[^q]*a$', list_colname)
check_colname('^(chl)[^q]*b$', list_colname)

test <- read_excel(i) %>% janitor::clean_names()

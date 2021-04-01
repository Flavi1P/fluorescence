require(sf)
require(tidyverse)
require(raster)
require(tidync)

ncpath <- "Data/BD6901481_001D.nc"
nc_path <- i
nc <- nc_open("Data/BD6901481_001D.nc")

VarName <- function(ncpath){
  t <- tidync(ncpath)
  tib_ <- as.tibble(t[[3]])
  myvar <- filter(tib_, grid == 'D10,D8') %>%
    pull(variables) %>%
    unlist()
  return(myvar)
}

variables <- VarName(nc_path)

get_qc <- function(varname, nc_path){
  nc <- nc_open(nc_path)
  t <- ncvar_get(nc, paste(varname, 'QC', sep = '_'))
  qc_vec <- unlist(str_split(t[[3]], pattern = ''))
  qc_vec <- gsub(' ', NA, qc_vec)
  nc_close(nc)
  return(qc_vec)
}

get_var <- function(varname, nc_path){
    var_tbl <- raster(nc_path,
                      varname = varname, na.rm = TRUE) %>%
      as.data.frame(xy = TRUE) %>%
      na.omit()
    
    names(var_tbl) <- c('depth', 'unuse', varname)
    return(var_tbl)
  
  
}

get_date <- function(nc){
  juld <- ncvar_get(nc, "JULD")
  juldqc <- ncvar_get(nc, "JULD_QC")
  origin<-NA
  origin<-as.POSIXct("1950-01-01 00:00:00", order="ymdhms") #convert juld->time
  time<-NA
  time<-origin + juld*3600*24
  time <- date(time)
  jd_qc<-NA
  jd_qc<-substr(ncvar_get(nc,"JULD_QC"),1,1)
  return(time)
  
}

nc <- nc_open(i)
chla <- ncvar_get(nc, 'CHLA')
names(nc$var)

extract_sd <- function(nc_path, vars){
  nc <- nc_open(nc_path)
  long_df <- data.frame('depth' = numeric(), 'variable' = character(), value = numeric())
  for(i in vars){
    var <- ncvar_get(nc, i)
    depth <- seq(1, length(var))
    table <- data.frame('depth' = depth, 'variable' = i, 'value' = var)
    long_df <- bind_rows(long_df, table)
  }
  lon <- ncvar_get(nc, 'LONGITUDE')
  lat <- ncvar_get(nc, 'LATITUDE')
  date <- get_date(nc)
  final_df <- long_df %>% pivot_wider(names_from = 'variable', values_from = 'value') %>% 
    mutate('date' = date,
           'lon' = lon,
           'lat' = lat)
  return(final_df)
}

extract_sd(i, c('CHLA', 'CHLA_ADJUSTED', 'CDOM'))


extract_bd <- function(nc_path, vars){
  variables <- VarName(nc_path)
  vars <- toupper(vars)
  long_df <- data.frame('depth' = numeric(), 'unuse' = numeric(), 'variable' = character(), value = numeric())
  for(i in vars){
    if(i %in% variables){
      temporary_df <- get_var(i, nc_path) %>% pivot_longer(3, names_to = 'variable', values_to = 'value')
      long_df <- bind_rows(long_df, temporary_df)
    }
    else{
      print(paste(i, 'is not a valid variable with this float choose from', list(tolower(variables)), sep = ' '))
      break()
    }
  }
  
  final_df <- pivot_wider(long_df, names_from = 'variable', values_from = 'value')
  return(final_df)
}

test <- extract_bgc(nc_path, c('chla', 'chla_adjusted', 'sal'))

ggplot(test) + 
  geom_path(aes(x = CDOM, y = -depth, colour = CDOM_QC))

read_hplc <- function(path){
  data <- read_excel(path) %>% janitor::clean_names()
  
 # cruise_name <- cruises_data[i,2]
  
  date <- select(data, matches('(^.+date)|(^date$)'))
  if(length(date) != 1){
    if(length(date) == 2){
      date <- select(date, matches('utc'))
    }
    if(length(date) == 0){
      year <- select(data, matches('year'))
      month <- select(data, matches('month'))
      day <- select(data, matches('day'))
      date <- paste(year[[1]], month[[1]], day[[1]], sep = '-')
      date <- data.frame("date" = date)
    }
  }
  
  lat <- select(data, matches('lat'))
  if(length(lat)!=1){
    lat <- lat[,1]
    if(length(lat) != 1){
      print(path)
      break()
    }
  }
  
  lon <- select(data, matches('lon'))
  if(length(lon) != 1){
    lon <- lon[,1]
    if(length(lon) != 1){
      print(path)
      break()
    }
  }
  
  ctd <- select(data, matches('ctd'))
  if(length(ctd) != 1){
    print(path)
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
  if(path != 14){
    
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
      chl_c3 <- NA
    }
    
    peri <- select(data, matches('peridinin'))
    if(length(peri) !=1){
      peri <- NA
    }
    
    but <- select(data, matches('(?=.*but)(?!.*qa)', perl = TRUE))
    if(length(but) !=1){
      but <- NA
    }
    
    fuco <- select(data, matches('^(?=fuco)(?!.*qa)', perl = TRUE))
    if(length(but) !=1){
      fuco <- NA
    } 
    
    neox <- select(data, matches('^(?=neox)(?!.*qa)', perl = TRUE))
    if(length(neox) !=1){
      neox <- NA
    }
    
    prasi <- select(data, matches('^(?=prasi)(?!.*qa)', perl = TRUE))
    if(length(prasi) !=1){
      prasi <- NA
    }
    
    viola <- select(data, matches('^(?=viola)(?!.*qa)', perl = TRUE))
    if(length(viola) !=1){
      viola <- NA
    }
    
    hex <- select(data, matches('(?=hex)(?!.*qa)', perl = TRUE))
    if(length(hex) !=1){
      hex <- NA
    }
    
    diad <- select(data, matches('(?=diad)(?!.*qa)', perl = TRUE))
    if(length(diad) !=1){
      diad <- NA
    }
    
    allo <- select(data, matches('(?=allo)(?!.*qa)', perl = TRUE))
    if(length(allo) !=1){
      allo <- NA
    }
    
    diat <- select(data, matches('(?=diat)(?!.*qa)', perl = TRUE))
    if(length(diat) !=1){
      diat <- NA
    }
    
    zea <- select(data, matches('(?=zea)(?!.*qa)', perl = TRUE))
    if(length(allo) !=1){
      allo <- NA
    }
    
    lutein <- select(data, matches('^(?=lutein)(?!.*qa)', perl = TRUE))
    if(length(lutein) !=1){
      lutein <- NA
    }
    
    dv_chlb <-  select(data, matches('(?=div.*b$)(?!.*qa)', perl = TRUE))
    if(length(dv_chlb) !=1){
        dv_chlb <- NA
      }
    
    chlb <- select(data, matches('^(?=chl.*b$)(?!.*qa)', perl = TRUE))
    if(length(chlb) !=1){
        chlb <- NA
    }
    tchlb <- select(data, matches('^(?=t.*chl.*b$)(?!.*qa)', perl = TRUE))
    if(length(tchlb) !=1){
      tchlb <- NA
    }
    
    dv_chla <- select(data, matches('(?=div.*a$)(?!.*qa)', perl = TRUE))
    if(length(dv_chla) !=1){
      dv_chla <- NA
    }
    
    chla <- select(data, matches('^(?!chlorophyllid)(?=chl.*a$)(?!.*qa)', perl = TRUE))
    if(length(chla) !=1){
     chla <- NA
    }
    
    tchla <- select(data, matches('^(?=t.*chl.*a$)(?!.*qa)', perl = TRUE))
    if(length(tchla) !=1){
      tchla <- NA
    }
    
    phaeophytin <- select(data, matches('(?=phaeophytin)(?!.*qa)', perl = TRUE))
    if(length(phaeophytin) !=1){
      phaeophytin <- NA
    }
    
  }
  
  
  data_interm <- tibble('date' = as.character(date[[1]]),
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
  return(data_interm)
}

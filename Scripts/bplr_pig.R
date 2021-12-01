library(tidyverse)
library(readxl)
source("Scripts/read_hplc.R")

argo <- read_csv("Output/Data/argo_matchup.csv")

bplr_argo <- filter(argo, code == "BPLR")

bplr_hplc <- read_hplc("Data/raw_excel/GreenEdge-Amundsen-pigments_flotteurs-300117.xlsx")

bplr_hplc_tomatch <- select(bplr_hplc, lon, lat, depth, peri, but, fuco, hex, allo, zea, chlb, t_chla, chla) %>% 
  mutate(lon = round(lon, 2),
         lat = round(lat, 2))

bplr_argo_tomatch <- select(bplr_argo, float_num, lon, lat, depth) %>% 
  mutate(lon = round(lon, 2),
         lat = round(lat, 2),
          depth = round(depth))

bplr_argo_tomatch$depth[bplr_argo_tomatch$depth == 4] <- 5

bplr_match <- left_join(bplr_argo_tomatch, bplr_hplc_tomatch)
bplr_match[bplr_match == "LOD"] <- "0"

bplr_match <- bplr_match %>% mutate(across(where(is.character), as.numeric))

bplr_argo_cols <- select(bplr_argo, -(chla:chlb)) %>% 
  mutate(lon = round(lon, 2),
         lat = round(lat, 2),
         depth = round(depth))

bplr_final <- left_join(bplr_argo_cols, bplr_match)

argo_nobplr <- filter(argo, code != "BPLR")

argo_final <- bind_rows(argo_nobplr, bplr_final)

write_csv(argo_final, "Output/Data/argo_matchup.csv")

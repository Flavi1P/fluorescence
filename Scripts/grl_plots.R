library(tidyverse)
library(sf)
library(s2)
library(readxl)
library(ggrepel)
library(gridExtra)
library(grid)
library(janitor)
library(Metrics)
library(cowplot)
library(ggtern)

library(patchwork)

path <- "Data/Longhurst"

argo <- read_csv("Data/argo_matchup.csv")

ref <- read_csv("Data/ref.csv")

longhurst_sf <- read_sf(dsn = path.expand(path), quiet = TRUE)

names(longhurst_sf) <- c("code", "region", "geometry")

pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(argo),
                                     function(i) {st_point(as.numeric(argo[i,c("lon", "lat") ]))}), list("crs" = 4326))) 

pnts_trans <- st_transform(pnts_sf, 4326)
longhurst_trans <- st_transform(longhurst_sf, 4326)  

sf_use_s2(FALSE)
argo$code <- apply(st_intersects(longhurst_trans, pnts_trans, sparse = FALSE), 2, 
                   function(col) { 
                     longhurst_trans[which(col), ]$code
                   })

sarc <- argo %>% filter(code == "SARC") %>% group_by(pres) %>% summarise_all(mean) %>% ungroup()
sarc$code <- "SARC"
sarc$lovbio <- "lovbio_sarc"

argo <- argo %>% filter(code != "SARC")
argo <- bind_rows(argo, sarc)

argo <- argo %>% mutate(code = case_when(code == "MEDI" & lon > 18 ~ "EMED",
                                                  code == "MEDI" & lon <= 18 ~ "WMED",
                                                  code != "EMED" & code != "WMED" ~ code))

med <- filter(argo)

Q <- quantile(med$ratio, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(med$ratio)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

argo <- filter(argo, code != "WMED" | (code == "WMED" & ratio < up & ratio > 0))

argo_new <- argo %>% mutate(fluo_biased = fluo *2, ratio = fluo_biased/t_chla)

region_argo <- argo_new  %>% group_by(code) %>% summarise_at(vars(ratio), c(mean, sd), na.rm = TRUE) %>% ungroup()
names(region_argo) <- c("code", "mean", "sd")

region_argo$sd <- ifelse(region_argo$sd > region_argo$mean, region_argo$mean, region_argo$sd)

codref <- read_excel("Data/Longhurst_Province_Summary.xls", 
                     range = "A17:B70", col_names = FALSE)
names(codref) <- c("code", "region")
region_argo <- left_join(region_argo, codref)


map_vec <- map_data("world")
g1 <- ggplot(filter(argo, code != 'ANTA'))+
  geom_point(aes(x = lon, y = lat, fill = code), size = 4, shape = 21)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec, fill = "gray28")+
  theme_bw(base_size = 16)+
  xlab("Longitude (°E)")+
  ylab("Latitude (°N)")+
  scale_fill_brewer(palette = 'Set1')+
  coord_quickmap()

region_argo <- region_argo[-7 ,]
region_argo$code <- as.character(region_argo$code)
g2 <- ggplot(region_argo)+
  geom_col(aes(x = reorder(code, mean), y = mean, fill = code))+
  geom_errorbar(aes(x = code, ymin = mean - sd, ymax = mean + sd))+
  xlab("Oceanic province")+
  ylab("Fluo/Chla ratio")+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size=1, linetype = "longdash", inherit.aes = F, width = 1)+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 20)

g1+g2





ggplot(filter(argo, code == "MEDI"))+
  geom_point(aes(x = fluo, y = -pres))+
  geom_point(aes(x = t_chla, y = - pres, colour = "chla"))+
  facet_wrap(.~(float), scales = "free")

g3 <- ggplot(argo)+
  coord_tern()+
  stat_interpolate_tern(geom = "polygon", formula = value~x+y,
                        method = lm, n = 50,
                        breaks = seq(0,7, by = 1),
                        aes(x = micro, y = nano, z = pico, value = ratio, fill =..level..), expand = 1)+
  geom_point(aes(x= micro, y = nano, z = pico, colour = code), size = 3)+
  scale_color_brewer(palette = "Set1", name = "Oceanic region")+
  theme_bw(base_size = 20)+
  guides(colour = FALSE)+
  weight_percent()+
  scale_fill_gradient(name = "Fluo/Chla ratio", low = "lightgrey", high = "gray45")
g3

g4 <- ggdraw()+
  draw_plot(g1,0,0.5,1,0.5)+
  draw_plot(g2,0,0,1,0.5)

grid.arrange(g4,g3, ncol = 2)  

#AFC####

library(vegan)


argo_afc <- select(argo, fluo, chla, ratio, all_of(pigments), code) %>% 
  mutate(fluo = log1p(fluo),
         tchla = log1p(chla),
         summ = rowSums(select(., fluo, chla, ratio))) %>% na.omit()

AFC <- rda(select(argo_afc, chla, ratio, fluo), scale = TRUE)

scores <- data.frame(scores(AFC, choices = c(1,2), display = "site"))
argo_afc <- bind_cols(argo_afc, scores)

pigscore <- data.frame(scores(AFC, choices = c(1,2), display = "species"))



fitscore <- envfit(AFC, select(argo_afc, all_of(pigments)))
fitarrow <- as.data.frame(fitscore$vectors$arrows)

argo_afc <- filter(argo_afc, code != 'ANTA')

ggplot(argo_afc)+
  geom_point(aes(x = PC1, y = PC2, colour = code))+
  geom_segment(aes(x = 0, xend = PC1 * 0.5, y = 0, yend = PC2 * 0.5), data = pigscore)+
  geom_text(aes(x = PC1 * 0.5, y = PC2 * 0.5, label = rownames(pigscore)), data = pigscore)+
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = fitarrow, colour = "#33a02c")+
  geom_text_repel(aes(x = PC1, y = PC2, label=rownames(fitarrow), fontface = 2), data = fitarrow)+
  scale_colour_brewer(palette = 'Set1')+
  theme_bw()+
  xlab('PC1')+
  ylab('PC2')




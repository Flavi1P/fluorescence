library(tidyverse)
library(janitor)
library(readxl)
library(patchwork)
library(vegan)
library(ggrepel)
library(treemap)
library(sf)
library(treemapify)
library(zoo)
library(caTools)
library(broom)
library(RColorBrewer)
library(lmodel2)

map_vec <- map_data("world")

types <- c('c', 'c', rep('n', 445))
lov_afc <- read_csv('Data/absorption/lov_afc.csv', col_types = as.list(types))

argo_hplc <- read_csv('Data/Raw HPLC Argo/hplc_argo_campaign') %>% rename("campagne" = cruise_name, "station" = station_name)

test <- bind_rows(lov_afc, argo_hplc)

ggplot(lov_afc)+
  geom_point(aes(x = lon, y = lat, colour = campagne))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  guides(colour = FALSE) +
  coord_quickmap()

lov_afc$asol440 = lov_afc$a440 + (0.0525 * (lov_afc$t_chla^0.855))
lov_afc$q_effect = lov_afc$x440 / lov_afc$asol440

lov_afc <- filter(lov_afc, q_effect <= 1)

ggplot(lov_afc)+
  geom_density(aes(x = q_effect))

lov_afc <- lov_afc %>% mutate(chla_470 = real470/t_chla,
                              chla_440 = real440/t_chla,
                              campagne = case_when(campagne %in% c("Pomme1 Leg 1", "Pomme 1 Leg 2") ~ "Pomme 1",
                                        campagne %in% c("Pomme2 Leg1", "Pomme 2 Leg 2") ~ "Pomme 2",
                                        campagne %in% c("Pomme3 Leg 1", "Pomme 3 Leg 2") ~ "Pomme 3",
                                        TRUE ~ campagne)) %>% 
  filter(z_zeu < 1)

#lov_afc <- select(lov_afc, campagne:real_440_470, x440, x470, a440, a470, real400:real550, pur440, pur470)

path = "Data/Longhurst"
longhurst_sf <- read_sf(dsn = path.expand(path), quiet = TRUE)

names(longhurst_sf) <- c("code", "region", "geometry")

#longhurst_sf %>% ggplot() + geom_sf(aes(fill = code))

pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(lov_afc),
                                     function(i) {st_point(as.numeric(lov_afc[i,c("lon", "lat") ]))}), list("crs" = 4326))) 
pnts_trans <- st_transform(pnts_sf, 4326)
longhurst_trans <- st_transform(longhurst_sf, 4326)  
sf_use_s2(FALSE)
lov_afc$code <- apply(st_intersects(longhurst_trans, pnts_trans, sparse = FALSE), 2, 
                      function(col) { 
                        longhurst_trans[which(col), ]$code
                      })
lov_afc$code <- as.character(lov_afc$code)


lov_map <- lov_afc %>% select(lon, lat, code, campagne) %>% distinct() %>% na.omit()

ggplot(lov_map)+
  geom_point(aes(x = lon, y = lat, colour = code))+
  geom_polygon(data = map_vec, aes(x = long, y = lat, group = group))+
  coord_quickmap()

# lov_afc <- filter(lov_afc, depth < 150 & t_chla < 1 & code != 'BENG' & code != 'character(0)' & code != 'CNRY') %>% 
#   mutate(code = case_when(code == "MEDI" & lon > 18 ~ "EMED",
#                           code == "MEDI" & lon <= 18 ~ "WMED",
#                           code != "EMED" & code != "WMED" ~ code))

x_axis <- seq(400, 520, by = 2)
weight440 <- dnorm(x_axis, 434, 8)
weight440 <- weight440 * 1/max(weight440)
weight470 <- dnorm(x_axis, 466, 11.5)
weight470 <- weight470 * 1/max(weight470)

lov_long <- lov_afc %>% mutate(row_id = c(1:nrow(lov_afc))) %>% select(row_id, x400:x520) %>%
  pivot_longer(x400:x520, names_to = 'wl', values_to = 'abs') %>% 
  group_by(row_id) %>% 
  mutate(weighted_440 = weighted.mean(x = abs, w = weight440),
         weighted_470 = weighted.mean(x = abs, w = weight470)) %>% 
  select(-abs, - wl) %>% 
  distinct()

lov_afc$zea[is.na(lov_afc$zea)] <- 0

lov_afc <- lov_afc %>% mutate(row_id = c(1:nrow(lov_afc))) %>% left_join(lov_long) %>%
  mutate(unquench_470 = weighted_470 * q_effect,
         unquench_440 = weighted_440 * q_effect,
         ratio = (zea + diadino)/(fuco + allo + x19bf + x19hf + peri + t_chlb + zea + diadino),
         ps_470 = unquench_470 - (unquench_470 * ratio),
         ps_440 = unquench_440 - (unquench_440 * ratio)) %>% 
  filter(code != "character(0)")

fitted_models <-  lov_afc %>% group_by(code) %>% do(model = lm(unquench_440 ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "t_chla") %>% 
  mutate(lambda = 440,
         slope = (estimate - mean(estimate)) / sd(estimate))

fitted_models_470 <-  lov_afc %>% group_by(code) %>% do(model = lm(weighted_470 ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "t_chla") %>% 
  mutate(lambda = 470,
         slope = (estimate - mean(estimate)) / sd(estimate))

stat_models_470 <-  lov_afc %>% group_by(code) %>% do(model = lm(weighted_470 ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, glance)) %>% 
  unnest(coef)

fitted_ps_470 <-  lov_afc %>% group_by(code) %>% filter(!is.na(ps_470)) %>% do(model = lm(ps_470 ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "t_chla") %>% 
  mutate(lambda = 470,
         slope = (estimate - mean(estimate)) / sd(estimate))

fitted_models_comb <-  lov_afc %>% group_by(code) %>% do(model = lm(weighted_470 + weighted_440 ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "t_chla") %>% 
  mutate(lambda = 470,
         slope = (estimate - mean(estimate)) / sd(estimate))

campagne_models_470 <-  lov_afc %>% group_by(campagne) %>% do(model = lm(weighted_470 ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(campagne, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "t_chla") %>% 
  mutate(lambda = 470,
         slope = (estimate - mean(estimate)) / sd(estimate))

models <- bind_rows(campagne_models_470, fitted_models_470) %>% select(slope, estimate, code, lambda) 

nb.cols <- length(unique(lov_afc$code))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
mycolors_campagne <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(lov_afc$campagne)))

ggplot(campagne_models_470)+
  geom_histogram(aes(x = reorder(campagne, estimate), y = estimate, fill = campagne), stat = 'identity')+
  geom_errorbar(aes(x = reorder(campagne, estimate), ymin = estimate - std.error, ymax = estimate + std.error))+
  ylab('slope')+
  xlab('oceanic province')+
  scale_fill_manual(values = mycolors_campagne)

ggplot(filter(fitted_models_470, code %in% c("ANTA", "ARCH", "ARCT", "BPLR", "MEDI", "SANT", "SPSG")))+
  geom_histogram(aes(x = reorder(code, estimate), y = estimate, fill = code), stat = 'identity')+
  geom_errorbar(aes(x = reorder(code, estimate), ymin = estimate - std.error, ymax = estimate + std.error))+
  ylab('slope')+
  xlab('oceanic province')+
  scale_fill_manual(values = mycolors)

ggplot(filter(fitted_models, code != "NASE" & code != "NADR" & code != "SSTC" & code != "CNRY"))+
  geom_histogram(aes(x = reorder(code, estimate), y = estimate, fill = code), stat = 'identity')+
  geom_errorbar(aes(x = reorder(code, estimate), ymin = estimate - std.error, ymax = estimate + std.error))+
  ylab('specific absorption')+
  xlab('oceanic province')+
  scale_fill_manual(values = mycolors)+
  theme_minimal()

gmap <- ggplot(lov_map)+
  geom_point(aes(x = lon, y = lat, fill = code), size = 3, shape = 21)+
  geom_polygon(data = map_vec, aes(x = long, y = lat, group = group))+
  coord_quickmap()+
  xlab("Longitude (°E)")+
  ylab("Latitude (°N)")+
  scale_fill_manual(values = mycolors)+
  theme_bw(base_size = 16)

yname <- expression(atop("a*(470)",~(m^2~"(mg"~"chla)"^{"-1"})))

fitted_models_470 <- fitted_models_470 %>% mutate(biome = case_when(code == "SANT" ~ "Polar",
                                                                    code == "BENG" ~ "Coastal",
                                                                    code == "CNRY" ~ "Coastal",
                                                                    code == "ANTA" ~ "Polar",
                                                                    code == "NATR" ~ "Equatorial",
                                                                    code == "MEDI" ~ "Temperate",
                                                                    code == "CHIL" ~ "Coastal",
                                                                    code == "ISSG" ~ "Equatorial",
                                                                    code == "SSTC" ~ "Temperate",
                                                                    code == "SPSG" ~ "Equatorial",
                                                                    code == "WARM" ~ "Equatorial",
                                                                    code == "PEQD" ~ "Equatorial",
                                                                    code == "NASE" ~ "Temperate",
                                                                    code == "NADR" ~ "Temperate"))

#Create a custom color scale
polarcolor <- brewer.pal(4,"PuBu")[c(3,4)]
medcolor <- brewer.pal(4, "Greens")
equatcolor <- brewer.pal(5, "OrRd")

myColors <- c(polarcolor, medcolor, equatcolor)

names(myColors) <- c("SANT", "ANTA", "MEDI", "SSTC", "NASE", "NADR", "NATR", "ISSG", "SPSG", "WARM", "PEQD")

fitted_models_470 <- filter(fitted_models_470, biome != "Coastal")

fitted_models_470$code <- factor(fitted_models_470$code, levels = c("SANT", "ANTA", "MEDI", "SSTC", "NASE", "NADR", "NATR", "ISSG", "SPSG", "WARM", "PEQD"))
gmap <- ggplot(filter(lov_map, code %in% fitted_models_470$code))+
  geom_polygon(data = map_vec, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = lon, y = lat, fill = code), size = 2, shape = 21)+
  coord_quickmap()+
  xlab("Longitude (°E)")+
  ylab("Latitude (°N)")+
  scale_fill_manual(values = myColors, name = "Code")+
  theme_bw(base_size = 11)


gmap/
  ggplot(fitted_models_470)+
  geom_histogram(aes(x = code, y = estimate, fill = code), stat = 'identity')+
  geom_errorbar(aes(x = code, ymin = estimate - std.error, ymax = estimate + std.error))+
  ylab(yname)+
  xlab('Bioregion')+
  scale_fill_manual(values = myColors, guide = "none")+
  theme_bw(base_size = 11)+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  plot_layout(guide = "collect", widths = 30)+
  geom_vline(xintercept = c(2.5, 6.5)) +
  annotate(geom = 'text', x = 1.5, y = 0.08, label = 'Polar') +
  annotate(geom = 'text', x = 4.5, y = 0.08, label = 'Temperate') +
  annotate(geom = 'text', x = 9, y = 0.08, label = 'Subtropical')+
  plot_layout(guide = "collect", height = c(1,1), nrow = 2)


ggsave("Output/Figures/fig_3_abs_newdim.jpg", width = 20, height = 13, unit = "cm", dpi = "print")

ggsave("Output/Figures/fig_3_abs.jpg", width = 20, height = 15, unit = "cm", dpi = "print")



gmap/
  ggplot(filter(fitted_ps_470, code != "NASE" & code != "NADR" & code != "SSTC" & code != "CNRY"))+
  geom_histogram(aes(x = reorder(code, estimate), y = estimate, fill = code), stat = 'identity')+
  geom_errorbar(aes(x = reorder(code, estimate), ymin = estimate - std.error, ymax = estimate + std.error))+
  ylab('Photosynthetic absorption')+
  xlab('Oceanic province')+
  scale_fill_manual(values = mycolors, guide = FALSE)+
  theme_bw(base_size = 18)


ggplot(filter(fitted_models_470, code != "NASE" & code != "NADR" & code != "SSTC" & code != "CNRY"))+
  geom_histogram(aes(x = reorder(code, estimate), y = q_effect, fill = code), stat = 'identity')+
  geom_errorbar(aes(x = reorder(code, estimate), ymin = q_effect - sd, ymax = q_effect + sd))+
  ylab('Specific absorption')+
  xlab('Oceanic province')+
  scale_fill_manual(values = mycolors, guide = FALSE)+
  theme_bw(base_size = 18)

  
model_440 <- lm(weighted_440 ~ t_chla, data = lov_afc)
summary(model_440)

model_470 <- lm(t_chla ~ weighted_470, data = lov_afc)
summary(model_470)

model_comb <- lm(t_chla ~ weighted_470 + weighted_440, data = lov_afc)
summary(model_comb)




max(fitted_models_470$estimate)

library(FactoMineR)

lov_afc <- left_join(lov_afc, campagne_models_470)

df_acp <- lov_afc %>% mutate(ratio = t_chla/weighted_470) %>% 
  select(estimate, t_chlb, peri, fuco, x19hf, x19bf, allo, zea, t_chla, campagne, weighted_470, estimate) %>% 
  mutate(t_chlb = t_chlb/t_chla,
         peri = peri/t_chla,
         fuco = fuco/t_chla,
         x19hf = x19hf/t_chla,
         x19bf, x19bf/t_chla,
         allo = allo/t_chla,
         zea = zea/t_chla)
test <- PCA(df_acp[,1:8], quanti.sup = 1,  scale.unit = TRUE)

coord <- data.frame(test$ind$coord)
var_coord <- bind_rows(data.frame(test$var$coord), data.frame(test$quanti.sup$coord))

result_acp <- bind_cols(df_acp, coord)

result_acp <- left_join(result_acp, fitted_models)
ggplot(filter(result_acp))+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = campagne, size = estimate))+
  scale_colour_manual(values = mycolors_campagne)+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2), data = var_coord)

ggplot(filter(result_acp))+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = estimate))+
  scale_colour_viridis_c()+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2), data = var_coord)

Q <- quantile(result_acp$ratio, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(result_acp$ratio)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

result_acp <- filter(result_acp, ratio < up & ratio > low)

ggplot(result_acp)+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = campagne))+
  scale_colour_manual(values = mycolors_campagne)+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1*3, yend = Dim.2*3), data = var_coord)+
  geom_text(aes(x = Dim.1*3.2, y = Dim.2*3.2, label = rownames(var_coord)), data = var_coord)+
  ylab("Dim 2 (19%)")+
  xlab("Dim 1 (29%)")

ggplot(result_acp)+
  geom_point(aes(x = Dim.1, y = Dim.3, colour = estimate))+
  scale_colour_distiller(palette = "Blues")+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1*3, yend = Dim.3*3), data = var_coord)+
  geom_text(aes(x = Dim.1*3.2, y = Dim.3*3.2, label = rownames(var_coord)), data = var_coord)+
  ylab("Dim 3 (16%)")+
  xlab("Dim 1 (29%)")

# ggplot(lov_afc)+
#   geom_point(aes(x = lon, y = lat, colour = code))+
#   geom_polygon(data = map_vec, aes(x = long, y = lat, group = group))+
#   coord_quickmap()+
#   scale_colour_manual(values = mycolors)+
#   theme_bw()
# 
# fitted_slope <-  lov_afc %>% group_by(code) %>% do(model = lm(log(t_chla) ~ log(a440), data = .)) %>% ungroup() %>% 
#   transmute(code, coef = map(model, tidy)) %>% 
#   unnest(coef) %>% 
#   mutate(lambda = 440)
# 
# fitted_slope2 <-  lov_afc %>% group_by(code) %>% do(model = lm(log(t_chla) ~ log(a470), data = .)) %>% ungroup() %>% 
#   transmute(code, coef = map(model, tidy)) %>% 
#   unnest(coef) %>% 
#   mutate(lambda = 470)
# 
# fitted_slope3 <-  lov_afc %>% group_by(code) %>% do(model = lm(log(t_chla) ~ log(real470) + log(real440), data = .)) %>% ungroup() %>% 
#   transmute(code, coef = map(model, tidy)) %>% 
#   unnest(coef) %>% 
#   mutate(lambda = 450)
# 
# slopes <- bind_rows(fitted_slope, fitted_slope2) %>% filter(term != '(Intercept)') %>%  select(estimate, code, lambda, ) %>% filter(code != 'character(0)') %>% janitor::clean_names()
# 
# ggplot(slopes)+
#   geom_point(aes(y = estimate, x = code, group = code), stat = "identity")+
#   theme_classic()+
#   ylab('Slope Chla~a_ps')+
#   xlab('Excitation wavelength (nm)')
# 
# 
# 
# rsq <- bind_rows(fitted_models, fitted_models2) %>% select(r.squared, code, lambda, ) 
# 
# ggplot(rsq)+
#   geom_boxplot(aes(y = r.squared, x = lambda, group = lambda))+
#   theme_classic()+
#   ylab('Slope Chla~a_ps')+
#   xlab('Excitation wavelength (nm)')
# 
# big_model <- lm(log(t_chla)~log(real440), data = lov_afc)
# lov_afc$resid440 <- residuals(big_model)
# lov_afc$fitted440 <- fitted(big_model)
# big_model2 <- lm(log(t_chla)~log(real470), data = lov_afc)
# lov_afc$resid470 <- residuals(big_model2)
# lov_afc$fitted470 <- fitted(big_model2)
# big_model3 <- lm(log(t_chla)~log(real440) + log(real470), data = lov_afc)
# lov_afc$fitted450 <- fitted(big_model3)
# lov_afc$resid450 <- residuals(big_model3)
# 
# 
# 
# lov_afc <- mutate(lov_afc, diff = abs(resid440) - abs(resid470))
# 
# 
# anova(big_model, big_model3)
# 
# ggplot(lov_afc)+
#   geom_point(aes(x = log(t_chla), y = fitted450))
# 
# 
# lov_prof <- lov_afc %>% group_by(lon, lat) %>% summarise_all(mean) %>% ungroup()
# 
# ggplot(lov_prof)+
#   geom_point(aes(x = lon, y = lat, colour = diff))+
#   geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
#   coord_quickmap()+
#   scale_color_viridis_c()
# 
# ggplot(lov_afc)+
#   geom_point(aes(x = t_chla, y = real470, colour = '470nm'))+
#   facet_wrap(.~code, scales = 'free_x')



test <- filter(longhurst_sf, code %in% unique(lov_afc$code))
test <- left_join(test, fitted_models_470)
test <- select(test, estimate, geometry)
plot(test)
plot(longhurst_sf)

lon <- unlist(test$geometry[1])[1:111]
lat <- unlist(test$geometry[1])[112:222]

ggplot()+
  geom_polygon(aes(x = lon, y = lat))

separated_coord <- test %>%
  mutate(lat = unlist(map(test$geometry,1)),
         long = unlist(map(test$geometry,2)))

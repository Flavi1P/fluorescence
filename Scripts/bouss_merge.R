library(tidyverse)
library(broom)
library(patchwork)

aphy <- read_csv("Data/Boussole/aphy_bouss.csv")
pig <- read_csv("Data/Boussole/bouss_pig.csv")

aphy$cruise <- as.numeric(aphy$cruise)
aphy$day <- as.numeric(aphy$day)
aphy <- mutate(aphy, sample = paste(cruise, ctd, depth))

pig <- mutate(pig, sample = paste(cruise_number, ctd, depth))

sample_aphy <- unique(aphy$sample)
sample_pig <- unique(pig$sample)

cruise_aphy <- unique(aphy$cruise)
cruise_pig <- unique(pig$cruise_number)

test <- which(!sample_aphy %in% sample_pig)
missing_cruise <- which(!cruise_aphy %in% cruise_pig)
cruise_aphy[missing_cruise]
missing_cruise2 <- which(!cruise_pig %in% cruise_aphy)
cruise_pig[missing_cruise2]

bouss_total <- left_join(pig, aphy, by = c("cruise_number" = "cruise", "depth", "day", "ctd"))

table(is.na(bouss_total$aphy))

bouss_wide <- bouss_total %>% pivot_wider( names_from = "lambda", values_from = "aphy")

bouss_wide <- filter(bouss_wide, depth < 105)
ggplot(bouss_wide)+
  geom_point(aes(x = date, y = -depth))

round_any <- function(x, z){
  if(x%%z > 0.5 * z){
    x <- x + (0.5 * z)
  }
  y <- floor(round(x)/z)
  y <- y * z
  return(y)
}

bouss_wide$depth2 <- 0
for(i in 1:length(bouss_wide$depth)){
  if(bouss_wide$depth[i] < 10){
    bouss_wide$depth2[i] <- round_any(bouss_wide$depth[i], 5)
  }
  else{
    bouss_wide$depth2[i] <- round_any(bouss_wide$depth[i], 10)
  }
  
}

ggplot(bouss_wide)+
  geom_point(aes(x = date, y = -depth2))
  

bouss_wide <- bouss_wide %>% mutate(ratio = (zea + diad)/(fuco + allo + but + hex + peri + tchlb + zea + diad)) %>% 
  filter(ratio >= 0)

ggplot(bouss_wide)+
  geom_density(aes(x = ratio))

boussole <- read_csv("Data/Boussole/boussole.csv") %>% select(date, press, depth, fluo, fluo_volt) %>% filter(depth < 100)

ggplot(boussole)+
  geom_point(aes(x = date, y = - depth))

boussole$depth2 <- 0
for(i in 1:length(boussole$depth)){
  if(boussole$depth[i] < 10){
    boussole$depth2[i] <- round_any(boussole$depth[i], 5)
  }
  else{
    boussole$depth2[i] <- round_any(boussole$depth[i], 10)
  }
  
}

ggplot(boussole)+
  geom_point(aes(x = date, y = - depth2))

boussole <- select(boussole, - depth)
test <- left_join(bouss_wide, boussole)

x_axis <- seq(400, 520, by = 1)
weight440 <- dnorm(x_axis, 434, 8)
weight440 <- weight440 * 1/max(weight440)
weight470 <- dnorm(x_axis, 466, 11.5)
weight470 <- weight470 * 1/max(weight470)

ggplot()+
  geom_line(aes(y = weight470, x = x_axis, colour = "470 nm excitation"))+
  geom_line(aes(y = weight440, x = x_axis, colour = "440 nm exciation"))+
  theme_bw(base_size = 18)+
  scale_color_brewer(palette = "Set1")+
  ylab("Normalized excitation")+
  xlab("Wavelength")

#ggsave("Output/paper_fig/excitation.png", height = 5, width = 10)

boussole_long <- test %>% janitor::clean_names() %>% mutate(row_id = c(1:nrow(test))) %>% select(row_id, x400:x520) %>%
  pivot_longer(x400:x520, names_to = 'wl', values_to = 'abs') %>% 
  group_by(row_id) %>% 
  mutate(weighted_440 = weighted.mean(x = abs, w = weight440),
         weighted_470 = weighted.mean(x = abs, w = weight470)) %>% 
  select(-abs, - wl) %>% 
  distinct()

bouss_fluo <- test %>% mutate(row_id = c(1:nrow(test))) %>%
  left_join(boussole_long) %>%
  select(date, depth2, fluo, fluo_volt, ratio, weighted_470, weighted_440, tchla.x, micro, nano, pico, peri:chla) %>% 
  mutate(ps_440 = weighted_440 - (weighted_440 * ratio),
         ps_470 = weighted_470 - (weighted_470 * ratio)) %>% 
  rename("tchla" = tchla.x) %>% 
  na.omit()

bouss <- bouss_fluo 

bouss$month <- lubridate::month(bouss$date)
bouss$monthabb <- month.abb[bouss$month]
bouss$year <- lubridate::year(bouss$date)

bouss <- bouss %>%
  mutate(
    season = case_when(
      month %in%  9:11 ~ "Fall",
      month %in%  c(12, 1, 2)  ~ "Winter",
      month %in%  3:5  ~ "Spring",
      month %in%  c(6, 7, 8) ~ "Summer")) %>% 
  filter(season != 'NA')

bouss$season <- factor(bouss$season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))

fitted_abs <-  bouss %>% group_by(depth2, season) %>% do(model = lm(weighted_470 ~ tchla, data = .)) %>% ungroup() %>% 
  transmute(depth2, season, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "tchla")

stat_abs <-  bouss %>% group_by(depth2, season) %>% do(model = lm(weighted_470 ~ tchla, data = .)) %>% ungroup() %>% 
  transmute(depth2, season, coef = map(model, glance)) %>% 
  unnest(coef)

fitted_yield <-  bouss %>% group_by(depth2, season) %>% do(model = lm(fluo_volt ~ weighted_440, data = .)) %>% ungroup() %>% 
  transmute(depth2, season, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "weighted_440")

stat_yield <-  bouss %>% group_by(depth2, season) %>% do(model = lm(fluo_volt ~ weighted_440, data = .)) %>% ungroup() %>% 
  transmute(depth2, season, coef = map(model, glance)) %>% 
  unnest(coef) 

fitted_slope <- bouss %>% group_by(depth2, season) %>% do(model = lm(fluo ~ tchla, data = .)) %>% ungroup() %>% 
  transmute(depth2, season, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "tchla")

stat_slope <- bouss %>% group_by(depth2, season) %>% do(model = lm(fluo ~ tchla, data = .)) %>% ungroup() %>% 
  transmute(depth2, season, coef = map(model, glance)) %>% 
  unnest(coef) 

fitted_ps <- bouss %>% group_by(depth2, season) %>% do(model = lm(ps_470 ~ tchla, data = .)) %>% ungroup() %>% 
  transmute(depth2, season, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "tchla")

yield_abs <- tibble("depth2" = fitted_abs$depth2,
                    "season" = fitted_abs$season,
                    "abs" = fitted_abs$estimate,
                    "abs_sd" = fitted_abs$std.error,
                    "yield" = fitted_yield$estimate,
                    "yield_sd" = fitted_yield$std.error,
                    "slope" = fitted_slope$estimate,
                    "slope_sd" = fitted_slope$std.error,
                    "ps" = fitted_ps$estimate,
                    "ps_sd" = fitted_ps$std.error) %>% 
  filter(depth2 <= 60)

#write_excel_csv(yield_abs, "output/bouss_table")
#yield_abs$depth2 <- factor(yield_abs$depth2, c("80", "70", "60", "50", "40", "30", "20", "10", "5"))

yname <- expression(atop("a*(470)"~(m^2~"(mg chla)"^{"-1"})))

ggplot(yield_abs)+
  geom_line(aes(y = abs, x = -depth2))+
  geom_point(aes(y = abs, x = -depth2))+
  geom_errorbar(aes(ymin = abs - abs_sd, ymax = abs + abs_sd, x = -depth2))+
  coord_flip()+
  facet_wrap(.~season, nrow = 1)+
  ylab(yname)+
  xlab("Depth (m)")+
  ylim(0,0.07)+
  theme_bw(base_size = 11)

ggsave("Output/Figures/specific_abs_bouss.png", width = 20, height = 10, units = "cm", dpi = 300)
#ggsave("Output/paper_fig/specific_abs_bouss_v2.png", width = 10, height = 5)

yname2 <- expression(atop("Chla specific photosynthetic absorption at 470 nm"~(m^-2~"(mg"%.%"chla)"^{"-1"})))


ggplot(yield_abs)+
  geom_bar(aes(y = yield, x = depth2, fill = depth2), stat = "identity")+
  geom_errorbar(aes(x = depth2, ymin = yield - yield_sd, ymax = yield + yield_sd))+
  scale_fill_brewer(palette = "Paired", direction = -1, guide = "none")+
  xlab("Depth (m)")+
  ylab("Quantum Yield of fluorescence (RFU.m²)")+
  coord_flip()+
  facet_wrap(.~season, nrow = 1)+
  theme_bw(base_size = 14)

yname3 <- expression(atop(phi~("RFU"~m^-1)))

yield_abs$yield[yield_abs$depth2 == 40 & yield_abs$season == "Spring"] <- 12.2

ggplot(yield_abs)+
  geom_line(aes(y = abs(yield), x = -depth2))+
  geom_point(aes(y = abs(yield), x = -depth2))+
  geom_errorbar(aes(ymin = yield - yield_sd, ymax = yield + yield_sd, x = -depth2))+
  coord_flip()+
  facet_wrap(.~season, nrow = 1)+
  ylab(yname)+
  xlab("Depth (m)")+
  ylab(yname3)+
  ylim(0,19)+
  theme_bw(base_size = 11)


ggsave("Output/Figures/yield_bouss.jpg",  width = 20, height = 10, units = "cm", dpi = 300)
#ggsave("Output/paper_fig/yield_abs_v2.png", width = 8, height = 4)

# compute package index ---------------------------------------------------

spectre <- read_excel("Data/absorption/Spectres_annick.xlsx") %>% janitor::clean_names() %>% 
  filter(lambda %in% c(430:450))


bouss_pig <- select(bouss, peri, but, fuco, hex, diad, allo, zea, tchlb, chla)

#change df to matrix 
bouss_mat <- as.matrix(bouss_pig)
bouss_mat <- t(bouss_mat)

#clean spectre and create a df with photosynthetical spectra
spectre[is.na(spectre)] <- 0

#create matrix
spectre_mat <- as.matrix(select(spectre, peri, x19_bf, fuco, x19_hf, diad, allox, zea, chl_b, chl_a))

#matrix multiplication to get a* * [pig1] + a* * [pig2] + ...
result <- spectre_mat %*%  bouss_mat

result <- t(result)

abs_calc <- rowSums(result)

bouss$a_calc <- abs_calc

bouss$package_index <- bouss$weighted_440 / bouss$a_calc
bouss$package_index <- bouss$package_index * 4
bouss <- bouss %>% mutate(weight440_unpackaged = weighted_440 * package_index,
                          weight470_unpackaged = weighted_470 * package_index)



ggplot(yield_abs)+
  geom_bar(aes(y = slope, x = depth2, fill = depth2), stat = "identity")+
  geom_errorbar(aes(x = depth2, ymin = slope - slope_sd, ymax = slope + slope_sd))+
  scale_fill_brewer(palette = "Paired", direction = -1, guide = FALSE)+
  xlab("Depth")+
  ylab("Slope factor (Volts.mg.Chl-1)")+
  coord_flip()+
  facet_wrap(.~season, nrow = 1)+
  theme_bw()


ggplot(yield_abs)+
  geom_point(aes(x = yield, y = abs, colour = depth2))

ggplot(yield_abs)+
  geom_point(aes(x = slope, y = yield))

yield_abs <- yield_abs %>% mutate(slope_calc = yield * abs)

ggplot(yield_abs)+
  geom_bar(aes(y = slope, x = -depth2), stat = "identity")+
  coord_flip()+
  facet_wrap(.~ season)

ggplot(yield_abs)+
  geom_point(aes(x = slope, y = slope_calc))

yield_abs$depth2 <- as.character(yield_abs$depth2)
yield_abs$depth2 <- as.numeric(yield_abs$depth2)


bouss <- left_join(bouss, yield_abs)

write_csv(file = "Output/Data/boussole_merged.csv", bouss)



    # bouss$row_id <- c(1:3451)
# pomme_long <- bouss %>% janitor::clean_names() %>% mutate(row_id = c(1:3451)) %>% select(row_id, x400:x500) %>%
#   pivot_longer(x400:x500, names_to = 'wl', values_to = 'x') %>% 
#   group_by(row_id) %>% 
#   mutate(weighted_440 = weighted.mean(x = x, w = weight440),
#          weighted_470 = weighted.mean(x = x, w = weight470)) %>% 
#   select(-x, - wl) %>% 
#   distinct()
# 
# bouss <- left_join(bouss, pomme_long)
# 
# bouss <- bouss %>%
#   mutate(
#     layer = case_when(
#       z_ze <=  0.5 ~ "surface",
#       z_ze <=  1 & z_ze > 0.5  ~ "middle",
#       z_ze <= 1.5 & z_ze > 1  ~ "depth"),
#     real470 = weighted_470 * ratio,
#     real440 = weighted_440 * ratio) %>% 
#   filter(layer != 'NA')
# 
# bouss$layer <- factor(bouss$layer, level = c('surface', 'middle', 'depth'))
# bouss$layer2 <- factor(bouss$depth2, level = c(5,10,20,30,40,50,60,70,80,100))
# ggplot(bouss)+
#   geom_point(aes(x = log(tchla.x), y = log(real470), colour = '470nm'))+
#   geom_point(aes(x = log(tchla.x), y = log(real440), colour = '440nm'))+
#   facet_grid(layer2 ~ season, scale = 'free')+
#   theme_bw()
# library(broom)
# fitted_models <-  bouss %>% select(season, layer2, tchla.x, weighted_440, weighted_470) %>% na.omit() %>% group_by(season, layer2) %>% do(tidy(lm(tchla.x ~ weighted_440, data = .))) %>% ungroup() %>% 
#   filter(term != '(Intercept)') %>%  
#   select(season, layer2, 'estimate_440' = estimate)
# 
# fitted_models2 <-  bouss %>% select(season, layer2, tchla.x, weighted_440, weighted_470) %>% na.omit() %>% group_by(season, layer2) %>% do(tidy(lm(tchla.x ~ weighted_470, data = .))) %>% ungroup() %>% 
#   filter(term != '(Intercept)') %>%  
#   select(season, layer2, 'estimate_470' = estimate)
# 
# bouss <- left_join(bouss, fitted_models) %>% left_join(fitted_models2)
# 
# coeff_long <- bouss %>% select(season, layer2, estimate_440, estimate_470) %>% 
#   pivot_longer(3:4, names_to = 'wl', values_to = 'slope') %>% distinct()
# 
# ggplot(filter(coeff_long, wl == 'estimate_470' & layer2 != "NA"))+
#   geom_tile(aes(x = season, y = layer2, fill = slope))+
#   geom_label(aes(x = season, y = layer2, label = round(slope, 1)))+
#   facet_grid(layer2~season, scales = 'free')+
#   scale_fill_viridis_c()
# 
# ggplot(coeff_long)+
#   geom_boxplot(aes(x = wl, y = slope, colour = wl))+
#   geom_jitter(aes(x = wl, y = slope, colour = wl), size = 1.5)+
#   scale_colour_brewer(palette = "Dark2")+
#   theme_bw()
# test <- bouss %>% filter(layer2 == "5")
# 
# ggplot(test)+
#   coord_tern()+
#   geom_point(aes(x= micro, y = nano, z = pico, colour = season), size = 3)+
#   scale_color_viridis_d()+
#   theme_bw()+
#   facet_wrap(.~season)
# 
# 
# ggplot(bouss)+
#   geom_point(aes(x = tchla.x, y = weighted_470))
# ggplot(bouss)+
#   coord_tern()+
#   geom_point(aes(x= micro, y = nano, z = pico, colour = estimate_470), size = 3)+
#   scale_color_viridis_c()+
#   theme_bw()+
#   facet_grid(layer ~ season)
# 
# pomme_pig <- bouss %>% select(season, layer,chla,  tchlb, peri, fuco, hex, but, allo, zea) %>% 
#   mutate(pigsum = tchlb + peri + fuco + hex + but + allo + zea) %>% 
#   pivot_longer(4:10, names_to = 'pig', values_to = 'concentration') %>% 
#   mutate(percent = concentration/pigsum)
# 
# ggplot(pomme_pig)+
#   geom_bar(aes(x = pig, y = percent, fill = pig), stat = 'identity', position = 'dodge')+
#   coord_polar(theta = "x", direction=1 )+
#   theme_bw()+
#   scale_fill_viridis_d()+
#   xlab('')+
#   facet_grid(layer~season)





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
map_vec <- read_csv("Data/map_vec")

types <- c('c', 'c', rep('n', 445))
lov_afc <- read_csv('Data/lov_afc.csv', col_types = as.list(types))

ggplot(lov_afc)+
  geom_point(aes(x = lon, y = lat, colour = campagne))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  guides(colour = FALSE) +
  coord_quickmap()

lov_afc <- lov_afc %>% mutate(chla_470 = real470/t_chla,
                              chla_440 = real440/t_chla)

rsq <- data.frame('chl' = seq(0.02, 1, 0.05),
                  'real_440' = NA,
                  'real_470' = NA,
                  'observed_440' = NA,
                  'observed_470' = NA,
                  'vitro_440' = NA,
                  'vitro_470' = NA)

for(i in c(1:length(rsq$chl))){
  threshold <- rsq$chl[i]
  df <- filter(lov_afc, t_chla < threshold)
  real_470 <- lm(t_chla~real470, data = df)
  real_440 <- lm(t_chla~real440, data = df)
  real_440 <- summary(real_440)$adj.r.squared
  real_470 <- summary(real_470)$adj.r.squared
  observed_470 <- summary(lm(t_chla~x470, data = df))$adj.r.squared
  observed_440 <- summary(lm(t_chla~x440, data = df))$adj.r.squared
  vitro_470 <- summary(lm(t_chla~a470, data = df))$adj.r.squared
  vitro_440 <- summary(lm(t_chla~a440, data = df))$adj.r.squared
  rsq$real_440[i] <- real_440
  rsq$real_470[i] <- real_470
  rsq$observed_440[i] <- observed_440
  rsq$observed_470[i] <- observed_470
  rsq$vitro_440[i] <- vitro_440
  rsq$vitro_470[i] <- vitro_470
}

rsq$diff <- rsq$real_440 - rsq$real_470
rsq$diffa <- rsq$vitro_440 - rsq$vitro_470


ggplot(rsq)+
  geom_path(aes(x = chl, y = diff, colour = 'aps'))+
  geom_point(aes(x = chl, y = diff))+
  geom_path(aes(x = chl, y = diffa, colour = 'a*ph'))+
  geom_point(aes(x = chl, y = diffa))+
  theme_bw()

ggplot(rsq)+
  geom_path(aes(x = chl, y = real_440, colour = 'Real 440'))+
  geom_path(aes(x = chl, y = real_470, colour = 'Real 470'))+
  geom_path(aes(x = chl, y = observed_440, colour = 'observed 440'))+
  geom_path(aes(x = chl, y = observed_470, colour = 'observed 470'))+
  geom_path(aes(x = chl, y = vitro_440, colour = 'vitro 440'))+
  geom_path(aes(x = chl, y = vitro_470, colour = 'vitro 470'))+
  scale_color_brewer(palette = 'Paired')+
  ylab('R² value')+
  xlab('Chla concentration')+
  ggtitle('R² of the relation aphy~chla')+
  xlim(0,1)

ggplot(rsq)+
  geom_path(aes(x = chl, y = vitro_440, colour = 'a*ph 440'), size = 2)+
  geom_path(aes(x = chl, y = vitro_470, colour = 'a*ph 470'), size = 2)+
  scale_color_brewer(palette = 'Set1')+
  ylab('R² value')+
  xlab('Chla concentration')+
  ggtitle('R² of the relation a*phy~chla')+
  xlim(0,1)+
  theme_bw(base_size = 20)

ggplot(rsq)+
  geom_path(aes(x = chl, y = real_440, colour = 'aps 440'), size = 2)+
  geom_path(aes(x = chl, y = real_470, colour = 'aps 470'), size = 2)+
  scale_color_brewer(palette = 'Set1')+
  ylab('R² value')+
  xlab('Chla concentration')+
  ggtitle('R² of the relation aps~chla')+
  xlim(0,1)+
  theme_bw(base_size = 20)

model <- lm(real440~t_chla, data = lov_afc)
test <- lov_afc %>% arrange(t_chla)

test <- lov_afc %>%
  filter(t_chla < 1) %>% 
  mutate(t_chla = plyr::round_any(t_chla, 0.1))
table(test$t_chla)

test <- test %>% group_by(t_chla) %>% 
  sample_n(size = 25)

test$cv440 <- runsd(test$real440, 40)/runmean(test$real440, 40)
test$cv470 <- runsd(test$real470, 40)/runmean(test$real470, 40)

test_summarised <- test %>% group_by(t_chla) %>% summarise_at(c('cv440', 'cv470'), mean)
ggplot(test_summarised)+
  geom_point(aes(x = t_chla ,y = cv440, colour = '440'))+
  geom_point(aes(x = t_chla, y = cv470, colour = '470'))

lov_model <- lov_afc %>% filter(t_chla < 3 & z_zeu <= 3)

ggplot(lov_model)+
  geom_point(aes(x = t_chla, y = real440), colour = 'skyblue2')+
  coord_trans(x = 'log', y = 'log')+
  ylab('aps 440')+
  theme_bw(base_size =20)+
  ggplot(lov_model)+
  geom_point(aes(x = t_chla, y = real470), colour = 'springgreen3')+
  coord_trans(x = 'log', y = 'log')+
  ylab('aps 470')+
  theme_bw(base_size = 20)

model_440 <- lm(t_chla~real440, data = lov_model)
model_470 <- lm(log(t_chla)~log(real470), data = lov_model)

summary(model_440)
summary(model_470)

AIC(model_440)
AIC(model_470)

exp((AIC(model_440)-AIC(model_470))/2)

df_model <- data.frame('chla' = lov_model$t_chla , 'resid440' = model_440$residuals, 'resid470' = model_470$residuals)

lov_afc$campagne_id <- sub("[1-9](.*)", "", lov_afc$campagne)

lov_afc %>% select(t_chla, real440, real470, campagne_id) %>% 
  pivot_longer(2:3, names_to = 'wl', values_to = 'abs') %>% 
  ggplot()+
  geom_point(aes(x = t_chla, y = abs, colour = wl))+
  coord_trans(x = 'log', y = 'log')+
  scale_colour_brewer(palette = 'Set1')+
  facet_wrap(.~campagne_id, scales = 'free_x')

lov_afc %>% select(t_chla, real440, real470, campagne_id) %>%
  filter(campagne_id == 'peacetime') %>% 
  pivot_longer(2:3, names_to = 'wl', values_to = 'abs') %>% 
  ggplot()+
  geom_point(aes(x = t_chla, y = abs, colour = wl))+
  coord_trans(x = 'log', y = 'log')+
  scale_colour_brewer(palette = 'Set1')


path = "Data/Longhurst"
longhurst_sf <- read_sf(dsn = path.expand(path), quiet = TRUE)

names(longhurst_sf) <- c("code", "region", "geometry")

#longhurst_sf %>% ggplot() + geom_sf(aes(fill = code))

pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(lov_afc),
                                     function(i) {st_point(as.numeric(lov_afc[i,c("lon", "lat") ]))}), list("crs" = 4326))) 
pnts_trans <- st_transform(pnts_sf, 4326)
longhurst_trans <- st_transform(longhurst_sf, 4326)  
lov_afc$code <- apply(st_intersects(longhurst_trans, pnts_trans, sparse = FALSE), 2, 
                   function(col) { 
                     longhurst_trans[which(col), ]$code
                   })
lov_afc$code <- as.character(lov_afc$code)
table(lov_afc$code)



rsq_camp <- data.frame('chl' = seq(0.02, 1, 0.05),
                            'biosope_440' = NA,
                            'biosope_470' = NA,
                         'peacetime_440' = NA,
                         'peacetime_470' = NA,
                         'soclim_440' = NA,
                         'soclim_470' = NA)


for(i in c(1:length(rsq_camp$chl))){
  threshold <- rsq$chl[i]
  dfbio <- filter(lov_afc, t_chla < threshold & campagne == 'Biosope')
  dfpeace <- filter(lov_afc, t_chla < threshold & campagne == 'peacetime')
  dfso <- filter(lov_afc, t_chla < threshold & campagne == 'soclim')
  biosope_470 <- summary(lm(t_chla~real470, data = dfbio))$adj.r.squared
  biosope_440 <- summary(lm(t_chla~real440, data = dfbio))$adj.r.squared
  peacetime_470 <- summary(lm(t_chla~real470, data = dfpeace))$adj.r.squared
  peacetime_440 <- summary(lm(t_chla~real440, data = dfpeace))$adj.r.squared
  soclim_470 <- summary(lm(t_chla~real470, data = dfso))$adj.r.squared
  soclim_440 <- summary(lm(t_chla~real440, data = dfso))$adj.r.squared
  
  rsq_camp$biosope_440[i] <- biosope_440
  rsq_camp$biosope_470[i] <- biosope_470
  
  rsq_camp$peacetime_440[i] <- peacetime_440
  rsq_camp$peacetime_470[i] <- peacetime_470
  
  rsq_camp$soclim_440[i] <- soclim_440
  rsq_camp$soclim_470[i] <- soclim_470
}

rsq_camp$diffbio <- rsq_camp$biosope_440 - rsq_camp$biosope_470
rsq_camp$diffpeace <- rsq_camp$peacetime_440 - rsq_camp$peacetime_470
rsq_camp$diffso <- rsq_camp$soclim_440 - rsq_camp$soclim_470


ggplot(rsq)+
  geom_path(aes(x = chl, y = diff, colour = 'aps'), size = 1.2)+
  geom_point(aes(x = chl, y = diff))+
  geom_path(aes(x = chl, y = diffa, colour = 'a*ph'), size = 1.2)+
  geom_point(aes(x = chl, y = diffa))+
  geom_path(aes(x = chl, y = diffbio, colour = 'aps biosope'), size = 1.2, data = rsq_camp)+
  geom_point(aes(x = chl, y = diffbio), data = rsq_camp)+
  geom_path(aes(x = chl, y = diffpeace, colour = 'aps peacetime'), size = 1.2, data = rsq_camp)+
  geom_point(aes(x = chl, y = diffpeace), data = rsq_camp)+
  geom_path(aes(x = chl, y = diffso, colour = 'aps soclim'), size = 1.2, data = rsq_camp)+
  geom_point(aes(x = chl, y = diffso), data = rsq_camp)+
  scale_color_brewer(palette = 'Set1')+
  theme_bw(base_size = 20)

lov_afc$seuil <- ceiling(lov_afc$chla/0.2)
rsq_threshold <- data.frame('chl' = seq(0.1, 1, 0.05),
                       'biosope_440' = NA,
                       'biosope_470' = NA,
                       'peacetime_440' = NA,
                       'peacetime_470' = NA,
                       'soclim_440' = NA,
                       'soclim_470' = NA)


for(i in c(1:length(rsq_threshold$chl))){
  threshold <- rsq_threshold$chl[i]
  dfbio <- filter(lov_afc, chla < threshold & campagne == 'Biosope')
  dfbio <- sample_n(dfbio, 25)
  dfpeace <- filter(lov_afc, chla < threshold & campagne == 'peacetime')
  dfpeace <- sample_n(dfpeace, 25)
  dfso <- filter(lov_afc, chla < threshold & campagne == 'soclim')
  dfso <- sample_n(dfso, 25)
  biosope_470 <- summary(lm(t_chla~real470, data = dfbio))$adj.r.squared
  biosope_440 <- summary(lm(t_chla~real440, data = dfbio))$adj.r.squared
  peacetime_470 <- summary(lm(t_chla~real470, data = dfpeace))$adj.r.squared
  peacetime_440 <- summary(lm(t_chla~real440, data = dfpeace))$adj.r.squared
  soclim_470 <- summary(lm(t_chla~real470, data = dfso))$adj.r.squared
  soclim_440 <- summary(lm(t_chla~real440, data = dfso))$adj.r.squared
  
  rsq_threshold$biosope_440[i] <- biosope_440
  rsq_threshold$biosope_470[i] <- biosope_470
  
  rsq_threshold$peacetime_440[i] <- peacetime_440
  rsq_threshold$peacetime_470[i] <- peacetime_470
  
  rsq_threshold$soclim_440[i] <- soclim_440
  rsq_threshold$soclim_470[i] <- soclim_470
}

rsq_threshold$diffbio <- rsq_threshold$biosope_440 - rsq_threshold$biosope_470
rsq_threshold$diffpeace <- rsq_threshold$peacetime_440 - rsq_threshold$peacetime_470
rsq_threshold$diffso <- rsq_threshold$soclim_440 - rsq_threshold$soclim_470

ggplot(rsq_threshold)+
  geom_path(aes(x = chl, y = diffbio, colour = 'aps biosope'), size = 1.2)+
  geom_point(aes(x = chl, y = diffbio))+
  geom_path(aes(x = chl, y = diffpeace, colour = 'aps peacetime'), size = 1.2)+
  geom_point(aes(x = chl, y = diffpeace))+
  geom_path(aes(x = chl, y = diffso, colour = 'aps soclim'), size = 1.2)+
  geom_point(aes(x = chl, y = diffso))+
  scale_color_brewer(palette = 'Set1')+
  theme_bw(base_size = 20)

rsq_camp2 <- rsq_camp %>% pivot_longer(2:7, names_to = 'wl', values_to = 'r') %>% 
  separate(wl, into = c('campagne', 'wl'), sep = '_')

ggplot(rsq_camp2)+
  geom_path(aes(x = chl, y = r, colour = wl), size = 2)+
  scale_color_brewer(palette = 'Set1')+
  ylab('R² value')+
  xlab('Chla concentration')+
  ggtitle('Adj R² aps~t_chla')+
  xlim(0,1)+
  theme_bw(base_size = 20)+
  facet_wrap(.~campagne, scales = 'free_y')

rsq_code <- data.frame('chl' = seq(0.02, 1, 0.05),
                       'medi_440' = NA,
                       'medi_470' = NA,
                       'nase_440' = NA,
                       'nase_470' = NA,
                       'spsg_440' = NA,
                       'spsg_470' = NA)


for(i in c(1:length(rsq_code$chl))){
  threshold <- rsq$chl[i]
  dfmedi <- filter(lov_afc, t_chla < threshold & code == 'MEDI')
  dfnase <- filter(lov_afc, t_chla < threshold & code == 'NASE')
  dfspsg <- filter(lov_afc, t_chla < threshold & code == 'SPSG')
  medi_470 <- summary(lm(t_chla~real470, data = dfmedi))$adj.r.squared
  medi_440 <- summary(lm(t_chla~real440, data = dfmedi))$adj.r.squared
  nase_470 <- summary(lm(t_chla~real470, data = dfnase))$adj.r.squared
  nase_440 <- summary(lm(t_chla~real440, data = dfnase))$adj.r.squared
  spsg_470 <- summary(lm(t_chla~real470, data = dfspsg))$adj.r.squared
  spsg_440 <- summary(lm(t_chla~real440, data = dfspsg))$adj.r.squared
  
  rsq_code$medi_440[i] <- medi_440
  rsq_code$medi_470[i] <- medi_470
  
  rsq_code$nase_440[i] <- nase_440
  rsq_code$nase_470[i] <- nase_470
  
  rsq_code$spsg_440[i] <- spsg_440
  rsq_code$spsg_470[i] <- spsg_470
}

rsq_code2 <- rsq_code %>% pivot_longer(2:7, names_to = 'wl', values_to = 'r') %>% 
  separate(wl, into = c('code', 'wl'), sep = '_')

ggplot(rsq_code2)+
  geom_path(aes(x = chl, y = r, colour = wl), size = 2)+
  scale_color_brewer(palette = 'Set1')+
  ylab('R² value')+
  xlab('Chla concentration')+
  ggtitle('Adj R² aps~t_chla')+
  xlim(0,1)+
  theme_bw(base_size = 20)+
  facet_wrap(.~code, scales = 'free_y')


ggplot(lov_afc)+
  geom_point(aes(x = lon, y = lat, colour = campagne_id))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()+
  theme_bw(base_size = 20)


lov_model <- lov_afc %>% filter(t_chla < 3 & z_zeu <= 3 & code %in% c('NASE', 'SPSG', 'MEDI'))

lov_model %>% select(t_chla, real440, real470, code, z_zeu) %>% 
  pivot_longer(2:3, names_to = 'wl', values_to = 'abs') %>% 
  ggplot()+
  geom_point(aes(x = t_chla, y = abs, colour = wl))+
  coord_trans(x = 'log', y = 'log')+
  facet_wrap(.~code)

rsq <- lov_model %>%
  mutate(layer = ceiling(z_zeu))%>%
  group_by(group) %>% 
  do(fit440 = lm(t_chla~real440, data = .),
     fit470 = lm(t_chla~real470, data = .))

rsq440 <- glance(rsq, fit440)
rsq470 <- glance(rsq, fit470)

rsq <- data.frame('group' = rsq440$group, 'r440' = rsq440$adj.r.squared, 'r470' = rsq470$adj.r.squared)



ggplot(rsq)+
  geom_tile(aes(x = group, y = layer, fill = r440 - r470))+
  scale_fill_distiller(palette = 'OrRd', direction = 1)+
  theme_minimal()+
  scale_x_discrete(position = 'top')

lov_model <- lov_afc %>% filter(t_chla < 3 & z_zeu <= 2)
lov_model[is.na(lov_model)] <- 0
diatom <- lov_model %>% mutate(pigsum = rowSums(select(., chla : tot_car)),
                                  diatom = zea/pigsum) %>% 
  group_by(campagne) %>% 
  summarise_at(vars(diatom), mean) %>% 
  ungroup()
which.max(diatom$diatom) #almofront-2


model470 <- lm(t_chla~real470, data = filter(lov_model, campagne == 'ALMOFRONT-2'))
model440 <- lm(t_chla~real440, data = filter(lov_model, campagne == 'ALMOFRONT-2'))

lov_innocent <- filter(lov_model, campagne %in% c('peacetime', 'Biosope', 'soclim')) %>%
  mutate(layer = ceiling(z_zeu))
pred470 <- predict(model470, lov_innocent)
pred440 <- predict(model440, lov_innocent)

rsq <- lov_innocent %>%
  mutate(layer = ceiling(z_zeu),
         pred440 = pred440,
         pred470 = pred470)%>%
  group_by(campagne, layer) %>% 
  do(fit440 = lm(t_chla~pred440, data = .),
     fit470 = lm(t_chla~pred470, data = .))

rsq440 <- glance(rsq, fit440)
rsq470 <- glance(rsq, fit470)

rsq <- data.frame('campagne' = rsq440$campagne, 'layer' = rsq440$layer, 'r440' = rsq440$adj.r.squared, 'r470' = rsq470$adj.r.squared)

ggplot(rsq)+
  geom_tile(aes(x = campagne, y = layer, fill = r440 - r470), colour = 'black')+
  scale_fill_distiller(palette = 'OrRd', direction = 1)+
  theme_minimal()+
  ylab('Z/Zeu')+
  scale_x_discrete(position = 'top')+
  scale_y_reverse()

lov_new <- lov_afc %>% filter(z_zeu < 3 & t_chla < 3)
model440 <- lm(real440~t_chla, data = lov_new)
model470 <- lm(real470~t_chla, data = lov_new)

lov_new$resid440 <- model440$residuals
lov_new$resid470 <- model470$residuals

lovbio <- lov_afc %>% filter(campagne == 'Biosope' & t_chla < 0.5) %>% 
  select(real440, real470, t_chla, lon, lat, depth) %>% 
  pivot_longer(c(1,2), names_to = 'wl', values_to = 'aps')


ggplot(lovbio)+
  geom_point(aes(x = t_chla, y = aps, colour = substr(wl, 5,7)))+
  geom_smooth(aes(x = t_chla, y = aps), method = 'lm', colour = 'black', se = FALSE)+
  scale_color_brewer(palette = 'Set1')+
  guides(colour = FALSE)+
  theme_bw(base_size = 20)+
  coord_trans(x = 'log', y = 'log')+
  facet_wrap(.~wl)
  


ggplot(lovbio)+
  geom_point(aes(x = lon, y = - depth, colour = resid470))+
  scale_color_viridis_c()


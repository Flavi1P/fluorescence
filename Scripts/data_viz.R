library(tidyverse)
library(FactoMineR)

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

argo <- read_csv("Data/merged_argo") %>%
  select(pigments, tchla, micro, nano, pico, ze, chla_adjusted) %>%
  mutate(campagne = "argo",
         fluo = chla_adjusted * 2,
         ratio = fluo/tchla) %>% 
  select(- chla_adjusted)
biosope <-  read_csv("Biosope/Data/biosope") %>%
  select(pigments, tchla, micro, nano, pico, ze, fluo_urel) %>%
  mutate(campagne = "biosope",
         fluo = fluo_urel,
         ratio = fluo/tchla) %>% 
  select(- fluo_urel)
boussole <- read_csv("Boussole/Data/boussole.csv") %>%
  select(pigments, tchla, micro, nano, pico, ze, fluo) %>%
  mutate(campagne = "boussole",
         ratio = fluo/tchla)


dataset_tall <- bind_rows(argo, biosope, boussole) %>% 
  select(micro, nano, pico, campagne, ratio) %>% 
  group_by(campagne) %>% 
  gather(key = "size_class", value = "frequence", 1:3) %>% 
  ungroup()


ggplot(dataset_tall)+
  geom_boxplot(aes(x = size_class, y = frequence, fill = campagne))+
  scale_fill_brewer(palette = "Set1", name = "Jeu de données", labels = c("argo" = "BGC-Argo", "biosope" = "BIOSOPE", "boussole" = "BOUSSOLE"))+
  ylim(0,1)+ylab("Proportion de la classe de taille")+ xlab("Classe de taille")+
  scale_x_discrete(labels = c("micro" = "Micro", "nano" = "Nano", "pico" = "Pico"))+
  theme_bw(base_size = 18)

#ggsave("Global_Dataset/Plots/size_class_boxplot.png", scale = 1.2)

ggplot(dataset_tall)+
  geom_violin(aes(x = campagne, y = ratio, fill = campagne))+
  scale_fill_brewer(palette = "Set1")+
  guides(fill = FALSE)+
  ylim(0,10)+
  ylab("Rapport [Chla]fluo / [Chla]HPLC")+
  xlab("Jeu de données")+
  scale_x_discrete(labels = c("argo" = "BGC-Argo", "biosope" = "BIOSOPE", "boussole" = "BOUSSOLE"))+
  theme_bw(base_size = 20)
#ggsave("Global_Dataset/Plots/ratio_violin.png")

dataset <- bind_rows(argo, biosope, boussole) %>%
  mutate(pigsum = rowSums(.[, 1:7])) %>%
  filter(pigsum != 0)


afc <- CA(select(dataset, pigments), row.sup = as.numeric(rownames(dataset[dataset$campagne != "argo",])))


scores <- data.frame(afc$row$coord)
scores_sup <- data.frame(afc$row.sup$coord)
scores <- bind_rows(scores, scores_sup)

dataset <- bind_cols(dataset,scores)
cos2 <- data.frame(afc$row$cos2)
cos2.sup <- data.frame(afc$row.sup$cos2)
cos2 <- bind_rows(cos2, cos2.sup)
dataset <- bind_cols(dataset, cos2)

pigscore <- data.frame(afc$col$coord)
dataset$cos <- (dataset$Dim.11+ dataset$Dim.21)

ggplot()+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = campagne), size = 0.9, data = dataset)+
  geom_segment(aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2), size = 0.8, data = pigscore)+
  ggrepel::geom_text_repel(aes(x = Dim.1, y = Dim.2, label = rownames(pigscore)), data = pigscore, size = 7)+
  scale_color_brewer(palette = "Set1", labels = c("argo" = "BGC-Argo", "biosope" = "BIOSOPE", "boussole" = "BOUSSOLE"), name = "Jeu de données")+
  ylab("CA2 (23%)") + xlab("CA1 (43%)")+
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)), alpha = guide_legend(title = "Cos²"))+
  theme(legend.position=c(.76,0.85),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.box = "horizontal",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                        colour = "grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "grey")) +
  coord_equal()


#ggsave("Global_Dataset/Plots/CA_comp.png", scale = 1.5)
biosope <- filter(dataset, campagne == "biosope") %>% 
  mutate(coord = paste(round(Dim.1,1), round(Dim.2,1), sep = ";"),
         dup = duplicated(coord)) %>% 
  filter(dup == "FALSE")

ggplot()+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = campagne), data = dataset)+
  geom_point(aes(x = Dim.1, y = Dim.2), data = biosope, colour = "blue")+
  geom_segment(aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2), data = pigscore)+
  geom_text(aes(x = Dim.1, y = Dim.2, label = rownames(pigscore)), data = pigscore)+
  scale_color_brewer(palette = "Set1")+
  theme_bw()

summary(lm(ratio~micro+nano+pico, data = filter(dataset, campagne == "argo")))
summary(lm(ratio~micro+nano+pico, data = filter(dataset, campagne == "biosope")))
summary(lm(ratio~micro+nano+pico, data = filter(dataset, campagne == "boussole" & ratio != "Inf")))

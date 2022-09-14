library(FactoMineR)
library(tidyverse)
library(ggrepel)

bouss <- read_csv("Output/Data/boussole_merged.csv")

bouss_acp <- bouss %>% mutate(tchlb = tchlb/tchla,
                              peri = peri/tchla,
                              fuco = fuco/tchla,
                              hex = hex/tchla,
                              but = but/tchla,
                              allo = allo/tchla,
                              zea = zea/tchla)


test <- PCA(bouss_acp[,c(8:19, 24)], quanti.sup = c(1,2,3,4), quali.sup = 13,  scale.unit = TRUE)

coord <- data.frame(test$ind$coord)
var_coord <- bind_rows(data.frame(test$var$coord))

supp_coord <- data.frame(test$quanti.sup$coord)
quali_coord <- data.frame(test$quali.sup$coord)

result_acp <- bind_cols(bouss, coord) %>% na.omit()

# Q <- quantile(result_acp$yield, probs=c(.25, .75), na.rm = FALSE)
# iqr <- IQR(result_acp$yield)
# 
# up <-  Q[2]+1.5*iqr # Upper Range
# low<- Q[1]-1.5*iqr # Lower Range
# 
# result_acp <- filter(result_acp, yield < up & yield > low)

# ggplot(result_acp)+
#   geom_point(aes(x = Dim.1, y = Dim.2, colour = abs))+
#   scale_colour_viridis_c()+
#   geom_segment(aes(x = 0, y = 0, xend = Dim.1*3, yend = Dim.2*3), data = var_coord)+
#   geom_text(aes(x = Dim.1*3.2, y = Dim.2*3.2, label = rownames(var_coord)), data = var_coord)+
#   ylab("Dim 2 (19%)")+
#   xlab("Dim 1 (29%)")



result_acp$depth2 <- factor(result_acp$depth2, c("5", "10", "20", "30", "40", "50", "60", "70", "80"))

result_acp <- filter(result_acp, depth2 %in% c("5", "10", "20", "30", "40", "50", "60"))

seasons_label <- result_acp %>% select(season, Dim.1, Dim.2) %>% 
  group_by(season) %>% 
  summarise_all(mean) %>%
  ungroup()

var_coord <- var_coord %>% mutate(pigname = rownames(.)) %>% 
  mutate(pigname = case_when(pigname == "but" ~ "19'-BF",
                             pigname == "hex" ~"19'-HF",
                             pigname == "fuco" ~"Fuco",
                             pigname == "allo" ~ "Allo",
                             pigname == "viola" ~"Viola",
                             pigname == "diad" ~ "Diad",
                             pigname == "zea" ~"Zea",
                             pigname ==  "peri"  ~"Peri"))


ggplot(result_acp)+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = depth2), size = 1.6)+
  scale_colour_brewer(palette = "Paired", name = "Depth (m)")+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1*2.7, yend = Dim.2*2.7), data = var_coord)+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1*3, yend = Dim.2*3), colour = "red", data = supp_coord)+
  geom_text_repel(aes(x = Dim.1*3, y = Dim.2*3, label = pigname), data = var_coord, size = 6)+
  geom_text_repel(aes(x = Dim.1*3.2, y = Dim.2*3.2, label = rownames(supp_coord)), colour = "red", data = supp_coord, size = 6)+
  geom_label(aes(x = Dim.1, y = Dim.2, label = season), data = seasons_label, size = 6)+
  ylab("PC2 (18.9%)")+
  xlab("PC1 (28.7%)")+
  xlim(-3, 3)+
  theme_bw(base_size = 11)

#ggsave("Output/paper_fig/acp_bouss.png", width = 9, height = 6)
#ggsave("Output/Figures/acp_bouss.jpg", dpi = "print", width = 20, height = 15, unit = "cm")

ggsave("Output/Figures/acp_bouss_redim.png", width = 9, height = 8)

colorname <- expression(atop("a*(470)",(m^2%.%"(mg"%.%"chla)"^{"-1"})))
astar <- expression(a[470]^{"*"})
aunit <- expression((m^-2%.%"(mg"%.%"chla)"^{"-1"}))

ggplot(result_acp)+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = abs))+
  scale_colour_viridis_c(name = colorname)+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1*2.7, yend = Dim.2*2.7), data = var_coord)+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1*3, yend = Dim.2*3), colour = "red", data = supp_coord)+
  geom_text_repel(aes(x = Dim.1*2.8, y = Dim.2*2.8, label = pigname), data = var_coord, size = 6)+
  geom_text_repel(aes(x = Dim.1*3.1, y = Dim.2*3.1, label = rownames(supp_coord)), colour = "red", data = supp_coord, size = 6)+
  ylab("PC2 (18.9%)")+
  xlab("PC1 (28.7%)")+
  xlim(-3, 3)+
  theme_bw(base_size = 11)

#ggsave("Output/paper_fig/acp_abs_bouss.png", width = 9, height = 6)
ggsave("Output/Figures/acp_abs_bouss.png", dpi = "print", width = 9, height = 8, units = "cm")

ggplot(result_acp)+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = slope))+
  scale_colour_viridis_c(name = "a*")+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1*3, yend = Dim.2*3), data = var_coord)+
  geom_segment(aes(x = 0, y = 0, xend = Dim.1*3, yend = Dim.2*3), colour = "red", data = supp_coord)+
  geom_text(aes(x = Dim.1*3.2, y = Dim.2*3.2, label = rownames(var_coord)), data = var_coord)+
  geom_text(aes(x = Dim.1*3.2, y = Dim.2*3.2, label = rownames(supp_coord)), colour = "red", data = supp_coord)+
  geom_label(aes(x = Dim.1, y = Dim.2, label = rownames(quali_coord)), data = quali_coord)+
  ylab("Dim 2 (18%)")+
  xlab("Dim 1 (22%)")+
  theme_bw()


# x_axis <- seq(400, 500, by = 1)
# weight440 <- dnorm(x_axis, 434, 8)
# weight440 <- weight440 * 1/max(weight440)
# weight470 <- dnorm(x_axis, 466, 11.5)
# weight470 <- weight470 * 1/max(weight470)
# 
# ggplot()+
#   geom_path(aes( x = x_axis, y = weight470, colour = "470"))+
#   geom_path(aes( x = x_axis, y = weight440, colour = "440"))+
#   theme_bw()
# 
# ggplot()+
#   geom_path(aes(x = x_axis, y = weight440, colour = 'abs 440'))+
#   geom_path(aes(x = x_axis, y = weight470, colour = 'abs 470'))+
#   scale_x_continuous(breaks = c(400, 420, 440, 460, 480, 500, 520))+
#   xlab('Wavelength (nm)')+
#   ylab('Weight')+
#   ggtitle('Approximation of excitation curves')+
#   theme_bw()
# 
# 
# bouss_ts_fit <- bouss %>% group_by(month) %>% do(model = lm(fluo ~ tchla, data = .)) %>% ungroup() %>% 
#   transmute(month, coef = map(model, tidy)) %>% 
#   unnest(coef) %>% 
#   filter(term == "tchla")
# 
# bouss_ts <- left_join(bouss, bouss_ts_fit)
# 
# bouss_ts$monthabb <- factor(bouss_ts$monthabb, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
# 
# ggplot(bouss_ts)+
#   geom_point(aes(x = monthabb, y = estimate))+
#   geom_line(aes(x = monthabb, y = estimate))
# 
# 
# ggplot(bouss_ts)+
#   geom_point(aes(x = date, y = estimate))+
#   geom_line(aes(x = date, y = estimate))+
#   ylim(0.15,0.5)
# 
# bouss_ts <- bouss %>% group_by(month, year, depth2) %>% summarise(chla_tot = mean(tchla), fluo_tot = mean(fluo))
# 
# bouss_ts <- left_join(bouss_ts, bouss)
# 
# bouss_ts$depth2 <- factor(bouss_ts$depth2, c("5", "10", "20", "30", "40", "50", "60", "70", "80"))
# 
# ggplot(bouss_ts)+
#   geom_line(aes(x = date, y = chla_tot, colour = depth2, group = depth2))+
#   facet_wrap(.~year, scale= "free_x", nrow = 3)+
#   scale_colour_brewer(palette = "Paired")+
#   theme_classic()

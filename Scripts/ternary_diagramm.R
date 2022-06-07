library(ggtern)
library(RColorBrewer)
argo <- read_csv("Output/Data/argo_matchup.csv")

#Create a custom color scale
polarcolor <- brewer.pal(5,"PuBu")[c(2,3,4,5)]
medcolor <- brewer.pal(4, "Greens")[c(3,4)]
equatcolor <- brewer.pal(3, "OrRd")[c(2,3)]

myColors <- c(polarcolor, medcolor, equatcolor)

names(myColors) <- c("SANT", "ARCT", "ANTA", "BPLR", "EMED", "WMED", "SPSG", "ARCH")

argo <- argo %>% mutate(wdp = 1.56 * fuco + 0.92 * peri + 4.8 * allo + 1.02 * but + 1.12 * hex + 1.51 * zea + 0.69 * chlb,
                        micro = (1.56 * fuco + 0.92 * peri)/wdp,
                        nano = (4.8 * allo + 1.02 * but + 1.51 * hex)/wdp,
                        pico = (1.51 * zea + 0.69 * chlb)/wdp)

argo <- filter(argo, lovbio != "lovbio024c" & lovbio != "lovbio026c" & lovbio != "lovbio027b" & lovbio != "lovbio040b")

argo <- argo %>% mutate(size_index = 0.1 * pico + 1 * nano + 10 * micro,
                        ratio = fluo_biased/t_chla) %>% 
  filter(ratio < 10 & size_index < 8)


ggtern(data = argo, aes(pico, nano, micro, value = ratio))+
  theme_rgbw(base_size = 18)+
  stat_interpolate_tern(geom="polygon",
                        formula=value~x+y,
                        method=lm,
                        aes(fill=..level..),
                        breaks = seq(1,7, by = 0.5),
                        expand=1,
                        )+
  geom_point(aes(colour = code))+
  scale_fill_gradient(low= "#ffffcc", high = "#41b6c4", name = "Slope factor")+
  scale_colour_manual(name = "code",values = myColors)

#ggsave("Output/paper_fig/ggtern.png", width = 10, height = 10)
ggsave("Output/Figures/ternary_diag.jpg", dpi = 300, width = 20, height = 15, unit = "cm")



ggplot(argo, aes(x = size_index, y = ratio))+
  geom_point(aes(colour = code))+
  geom_smooth(method = "lm")+
  scale_colour_brewer(palette = "Set1")
